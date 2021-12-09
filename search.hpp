#ifndef SEARCH_HPP
#define SEARCH_HPP

#include <functional>
#include <assert.h>


namespace search {
double binary(double start, double end, const double target, const double precision, const std::function<double(double)> fun)
{
	if (fun(start) - target > 0) {
		// switch start and end
		start += end;
		end = start - end;
		start -= end;
	}
	assert(fun(end) > 0);
	double midval, error, midpoint;
	do {
		midpoint = (start + end) / 2;
		if (midpoint == start || midpoint == end) {
			throw std::runtime_error("reached maximum precision");
			break; // maximum precision has been reached
		}
		midval = fun(midpoint);
		if (midval < 0)
			start = midpoint;
		else
			end = midpoint;
		error = std::abs(end - start); //std::abs(midval - target);
		//std::cout << "midpoint " << midpoint << ", error " << error << std::endl;
	} while (error > precision);
	return midpoint;
}
}
#endif