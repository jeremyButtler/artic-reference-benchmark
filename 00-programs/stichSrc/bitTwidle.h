/*#########################################################
# Name: bitTwidle
# Use:
#  - Holds bit twidling functions
# Libraries:
# C Standard Libraries:
# Sources:
#  - https://graphics.stanford.edu/~seander/bithacks.html
#########################################################*/

#ifndef BITTWIDLE_H
#define BITTWIDLE_H

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' bitTwidle SOH: Start Of Header
' Fun-01 bitMin:
'  - Uses bit twidling to do a branchless min. From my
'    experince this is a bit faster than an if min.
' Fun-02 bitMax:
'  - Uses bit twidling to do a branchless max. From my
'    experince this is a bit faster than an if max.
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Name: bitMin (Fun-01:)
| Use:
|  - Uses bit twidling to do a branchless min. From my
|    experince this is a bit faster than an if min.
| Input:
|  - ret:
|    o value to store min in. This can be x or y
|  - x:
|    o First value to check
|  - y:
|    o Second value to check. This is kept if X == Y
| Output:
|  - Modifies:
|    - ret to hold the minimum value or Y if both are equal
\--------------------------------------------------------*/
#define bitMin(ret, x, y) { \
   (ret) = (x) ^ ( ((x) ^ (y)) & ( -((x) > (y)) ) ); \
}

/*--------------------------------------------------------\
| Name: bitMax (Fun-02:)
| Use:
|  - Uses bit twidling to do a branchless max. From my
|    experince this is a bit faster than an if max.
| Input:
|  - ret:
|    o value to store max in. This can be x or y
|  - x:
|    o First value to check
|  - y:
|    o Second value to check. This is kept if X == Y
| Output:
|  - Modifies:
|    - ret to hold the maximum value or y if both are equal
\--------------------------------------------------------*/

#define bitMax(ret, x, y) { \
   (ret) = (x) ^ ( ((x) ^ (y)) & ( -((x) < (y)) ) ); \
}
#endif
