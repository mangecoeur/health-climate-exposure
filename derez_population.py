import rasterio
from affine import Affine

from config import POP_DATA_SRC

population_grid_2015_file = (POP_DATA_SRC / 'nasa_grid' /
                             'gpw-v4-population-count-adjusted-to-2015-unwpp-country-totals-2015' /
                             'gpw-v4-population-count-adjusted-to-2015-unwpp-country-totals_2015.tif')

with rasterio.open(str(population_grid_2015_file)) as pop:
    print(pop.meta)
    pop_meta = pop.meta
    trns = pop.transform
    population = pop.read(1, masked=True)

population.fill_value = 0
population = population.filled()

# Sum every other row
population = population[::2, :] + population[1::2, :]
# Sum every other column
population = population[:, ::2] + population[:, 1::2]

# Output affine scaled by 2
newtrans = Affine(trns.a * 2, trns.b, trns.c, trns.d, trns.e * 2, trns.f)

# NOTE: loop the above to reduce the resolution further.
# Reduction to 1/4 of the original size already makes life much easier

print('Saving')
print(population.shape)
with rasterio.open('/Users/jonathanchambers/pop_lowres.tif', 'w', driver='GTiff',
                   height=population.shape[0],
                   width=population.shape[1],
                   count=1,
                   dtype=population.dtype,
                   crs=pop_meta['crs'],
                   transform=newtrans,
                   compress='lzw') as new_dataset:
    new_dataset.write(population, 1)
# new_dataset.write_mask(np.invert(population.mask))
