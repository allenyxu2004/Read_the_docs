#include "cordiccart2pol.h"
#include <math.h>

data_t Kvalues[NO_ITER] = {1,	0.500000000000000,	0.250000000000000,	0.125000000000000,	0.0625000000000000,	0.0312500000000000,	0.0156250000000000,	0.00781250000000000,	0.00390625000000000,	0.00195312500000000,	0.000976562500000000,	0.000488281250000000,	0.000244140625000000,	0.000122070312500000,	6.10351562500000e-05,	3.05175781250000e-05};

data_t angles[NO_ITER] = {0.785398163397448,	0.463647609000806,	0.244978663126864,	0.124354994546761,	0.0624188099959574,	0.0312398334302683,	0.0156237286204768,	0.00781234106010111,	0.00390623013196697,	0.00195312251647882,	0.000976562189559320,	0.000488281211194898,	0.000244140620149362,	0.000122070311893670,	6.10351561742088e-05,	3.05175781155261e-05};

data_t PI = 3.141592653589793;

void cordiccart2pol(data_t x, data_t y, data_t * r,  data_t * theta)
{	
	// If in Quandrant II or III, reflect across the origin to land in Quandrant I or IV
	cos_sin_type current_x = x;
	cos_sin_type current_y = y;

	if (x < 0) {
		current_x = -x;
		current_y = -y;
		*theta = (y >= 0) ? PI : -PI;
	}
	else {
		*theta = 0;
	}
	
	cos_sin_type cordic_factor = 1;

	for (int i = 0; i < NO_ITER; i++) {
		// Determine if we want to rotate by a positive or negative angle
		int sigma = (current_y > 0) ? -1 : 1;

		cos_sin_type x_shift = current_x * sigma * Kvalues[i];
		cos_sin_type y_shift = current_y * sigma * Kvalues[i];

		current_x = current_x - y_shift;
		current_y = current_y + x_shift;

		*theta -= (sigma * angles[i]);

		// For now, directly calculate the cordic gain
		cordic_factor *= (sqrtf(1 + (Kvalues[i] * Kvalues[i])));
	}

	// Resulting x-component is the magnitude multiplied by the cordic gain
	*r = current_x / cordic_factor;
	
}
