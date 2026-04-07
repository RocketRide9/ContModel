#define real double

kernel void write_ids_1d(
	global int* ids,
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

kernel void fill(
	global real* dst,
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
		dst[i] = item;
	}
}
