#define real double
#if !defined(CUDA)
	#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#endif

__kernel void write_ids_1d(
	__global int* ids,
	int n
) {
	int i = get_global_id(0);
	if (i < 0) {
		printf("this is bad");
		return;
	}
	if (i < n)
	{
		ids[i] = i;
	}
}

__kernel void fill(
	__global real* dst,
	real item,
	int n
) {
	int i = get_global_id(0);
	if (i < 0) {
		printf("this is bad");
		return;
	}
	if (i < n)
	{
		dst[i] = item * 1.0;
	}
}

// y *= a
__kernel
void scale
(
	__global real* y,
	const real a,
	const int n
)
{
	uint i = get_global_id(0);
	if (i < n)
	{
		//real v = a;
		//y[i] += 0.5;
		y[i] *= 0.5;
	}
}