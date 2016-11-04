R"(

/**
 * @brief Source code for the kernel fuctions used by the program.
 *
 * @file
 * @include kernels.cl
 */

/**
 * @brief Multiplies two complex numbers.
 */
inline float2 complexMultiply(float2 a, float2 b)
{
	return (float2)(a.x*b.x - a.y*b.y, a.x*b.y + a.y*b.x);
}

/**
 * @brief Applies the filter.
 *
 * Multiplies the complex numbers in a by corresponding elements in b.
 */
__kernel void filter(__global float2* a, __global float2* b)
{
	int id0 = get_global_id(0);
	int id1 = get_global_id(1);
	int size1 = get_global_size(1);

	int id = id0*size1 + id1;

	a[id] = complexMultiply(a[id], b[id1]);
}

/**
 * @brief Assigns zero to all elements.
 */
__kernel void zero(__global float* a)
{
	int id0 = get_global_id(0);
	a[id0] = 0;
}

)"
