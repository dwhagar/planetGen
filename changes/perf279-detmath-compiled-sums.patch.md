**Headline:** the long sums of the seasonal temperature solve (`detmath.matmul`) are added in one compiled pass instead of arrays of products, with the same bits.

### Changed
- `detmath.matmul` adds a long sum (from 100,000 products) with a small numba kernel (`util/detmath_kernel.py`) that does what the loop and `numpy.add.accumulate` did: each output starts from the product at index 0 and adds the products at index 1, 2 and so on in order, each product and each addition correctly rounded once, no fused multiply-add and no reordering. The chains of the outputs advance together, so a sum runs at the speed of the adders and not of their latency. The seasonal solve's three sums per pass (a Gram row, a right-hand side and the mixed series, 6 series of up to 190,080 values) take 4 ms instead of 17 ms.
- Checked against the real build: 413 calls of `temperature_series` made while the rich test galaxy was built (pre-placed systems and promoted belt members) give the same temperatures, pass counts and convergence to the last bit with and without the kernel, and take 37.6 s instead of 65.8 s. Where numba cannot be imported or cannot cache, `matmul` keeps the numpy way and its bits; the kernel is imported inside `matmul`, so the web app still never imports numba.

generation-output: unchanged
