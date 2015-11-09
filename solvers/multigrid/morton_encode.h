#ifndef MORTON_ENCODE_H
#define MORTON_ENCODE_H

// SOURCE: taken from
// http://and-what-happened.blogspot.com/2011/08/fast-2d-and-3d-hilbert-curves-and.html

// TODO: make this uint instead?
int morton_encode(int i, int j, int k)
{
    // pack 3 10-bit to 32 bit value
    i &= 0x000003ff;
    j &= 0x000003ff;
    k &= 0x000003ff;

    i |= ( i << 16 );
    j |= ( j << 16 );
    k |= ( k << 16 );

    i &= 0x030000ff;
    j &= 0x030000ff;
    k &= 0x030000ff;

    i |= ( i << 8 );
    j |= ( j << 8 );
    k |= ( k << 8 );

    i &= 0x0300f00f;
    j &= 0x0300f00f;
    k &= 0x0300f00f;

    i |= ( i << 4 );
    j |= ( j << 4 );
    k |= ( k << 4 );

    i &= 0x030c30c3;
    j &= 0x030c30c3;
    k &= 0x030c30c3;

    i |= ( i << 2 );
    j |= ( j << 2 );
    k |= ( k << 2 );

    i &= 0x09249249;
    j &= 0x09249249;
    k &= 0x09249249;

    // make it contiguous in k, so kji
    return ( k | (j << 1) | (i << 2) );
}

void morton_decode(int morton, int *i, int *j, int *k)
{
    int val1 = morton;
    int val2 = ( val1 >> 1);
    int val3 = ( val1 >> 2);

    val1 &= 0x09249249;
    val2 &= 0x09249249;
    val3 &= 0x09249249;

    val1 |= (val1 >> 2);
    val2 |= (val2 >> 2);
    val3 |= (val3 >> 2);

    val1 &= 0x030c30c3;
    val2 &= 0x030c30c3;
    val3 &= 0x030c30c3;

    val1 |= (val1 >> 4);
    val2 |= (val2 >> 4);
    val3 |= (val3 >> 4);

    val1 &= 0x0300f00f;
    val2 &= 0x0300f00f;
    val3 &= 0x0300f00f;

    val1 |= (val1 >> 8);
    val2 |= (val2 >> 8);
    val3 |= (val3 >> 8);

    val1 &= 0x030000ff;
    val2 &= 0x030000ff;
    val3 &= 0x030000ff;

    val1 |= (val1 >> 16);
    val2 |= (val2 >> 16);
    val3 |= (val3 >> 16);

    val1 &= 0x000003ff;
    val2 &= 0x000003ff;
    val3 &= 0x000003ff;

    *i = val3;
    *j = val2;
    *k = val1;
}


#endif
