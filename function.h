/* Tumour cord growth model */

/*
 * Copyright (C) 2006 Sergey Astanin, Luigi Preziosi, Andrea Tosin
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#ifndef FUNCTION_H
#define FUNCTION_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <exception>
#include <vector>
#include <string>
#include <algorithm>
#include <iostream>

using std::exception;
using std::vector;
using std::string;
using std::sort;

class FunctionException : exception {
private:
	string msg;
public:
	~FunctionException() throw() {}
	FunctionException(string _msg) : msg(_msg) {}
	const char* what() const throw() {
		return msg.c_str();
	}
};

/** @brief ad-hoc constructable function object */
class ADoubleFunction {
public:
	virtual ~ADoubleFunction() {};
	virtual double operator()(double const arg) const = 0;
	virtual ADoubleFunction* clone() const = 0; // uses copy constructor
	double eval(double const arg) const { return this->operator()(arg); }
};

/** @brief sum of two functions */
class SumFunction : public ADoubleFunction {
private:
	ADoubleFunction* f1;
	ADoubleFunction* f2;
public:
	virtual ~SumFunction() {
		if (f1) {
			delete f1;
			f1=(ADoubleFunction*)0;
		}
		if (f2) {
			delete f2;
			f2=(ADoubleFunction*)0;
		}
	}

	SumFunction(const ADoubleFunction& _f1, const ADoubleFunction& _f2)
		: f1(_f1.clone()), f2(_f2.clone()) {}

	SumFunction(const ADoubleFunction* _pf1, const ADoubleFunction* _pf2)
	throw(FunctionException)
		: f1(_pf1->clone()), f2(_pf2->clone()) {
		if ((!_pf1)||(!_pf2)) {
			throw FunctionException("SumFunction::SumFunction"
			"(const ADoubleFunction*, const ADoubleFunction*): "
			"cannot build sum from null pointers to functions");
		}
	}

	SumFunction(const SumFunction& f)
		: f1(f.f1->clone()), f2(f.f2->clone()) {}
	
	SumFunction& operator=(const SumFunction& rhs) {
		if (this != &rhs) {
			ADoubleFunction* f=rhs.f1->clone();
			if (f1) {
				delete f1;
				f1=(ADoubleFunction*)0;
			}
			f1=f;
			f=rhs.f2->clone();
			if (f2) {
				delete f2;
				f2=(ADoubleFunction*)0;
			}
			f2=f;
		}
		return *this;
	}

	virtual double operator()(const double arg) const throw(FunctionException) {
		if ((!f1) || (!f2)) {
			throw FunctionException("SumFunction::operator(): "
				"null pointer");
		}
		return ((*f1)(arg)+(*f2)(arg));
	}

	virtual ADoubleFunction* clone() const {
		return (ADoubleFunction*)(new SumFunction(*this));
	}
};

/** @brief product of two functions */
class ProductFunction : public ADoubleFunction {
private:
	ADoubleFunction* f1;
	ADoubleFunction* f2;
public:
	virtual ~ProductFunction() {
		if (f1) {
			delete f1;
			f1=(ADoubleFunction*)0;
		}
		if (f2) {
			delete f2;
			f2=(ADoubleFunction*)0;
		}
	}

	ProductFunction(const ADoubleFunction& _f1, const ADoubleFunction& _f2)
		: f1(_f1.clone()), f2(_f2.clone()) {}

	ProductFunction(const ADoubleFunction* _pf1, const ADoubleFunction* _pf2)
	throw(FunctionException)
		: f1(_pf1->clone()), f2(_pf2->clone()) {
		if ((!_pf1)||(!_pf2)) {
			throw FunctionException("ProductFunction::"
			"ProductFunction(const ADoubleFunction*, "
			"const ADoubleFunction*): "
			"cannot build sum from null pointers to functions");
		}
	}

	ProductFunction(const ProductFunction& f)
		: f1(f.f1->clone()), f2(f.f2->clone()) {}
	
	ProductFunction& operator=(const ProductFunction& rhs) {
		if (this != &rhs) {
			ADoubleFunction* f=rhs.f1->clone();
			if (f1) {
				delete f1;
				f1=(ADoubleFunction*)0;
			}
			f1=f;
			f=rhs.f2->clone();
			if (f2) {
				delete f2;
				f2=(ADoubleFunction*)0;
			}
			f2=f;
		}
		return *this;
	}

	virtual double operator()(const double arg) const throw(FunctionException) {
		if ((!f1) || (!f2)) {
			throw FunctionException("ProductFunction::operator(): "
				"null pointer");
		}
		return ((*f1)(arg)*(*f2)(arg));
	}

	virtual ADoubleFunction* clone() const {
		return (ADoubleFunction*)(new ProductFunction(*this));
	}
};

class ConstFunction : public ADoubleFunction {
private:
	double value; ///< constant value
public:
	virtual ~ConstFunction() {}
	ConstFunction(const double _value) : value(_value) {}
	virtual double operator()(const double arg) const { return value; }
	virtual ADoubleFunction* clone() const {
		return (ADoubleFunction*)(new ConstFunction(*this));
	}
};

class LinearFunction : public ADoubleFunction {
private:
	double slope; ///< slope tan
	double yvalue; ///< intersection with y-axis
public:
	virtual ~LinearFunction() {}
	LinearFunction(const double _slope=1.0, const double _yvalue=0.0)
		: slope(_slope), yvalue(_yvalue) {}
	virtual double operator()(const double arg) const { return arg*slope+yvalue; }
	virtual ADoubleFunction* clone() const {
		return (ADoubleFunction*)(new LinearFunction(*this));
	}
};

class SigmaFunction : public ADoubleFunction {
private:
	double arg0; ///< Sigma's zero
	double kright; ///< right slope
	double kleft; ///< left slope

public:
	virtual ~SigmaFunction() {};

	SigmaFunction(double const phi0, double const k1, double const s1)
		: arg0(phi0), kright(k1), kleft(s1) {}

	virtual double operator()(double const arg) const {
		if (arg > arg0) {
			return kright*(arg-arg0);
		} else {
			return kleft*(arg-arg0);
		}
	}

	virtual ADoubleFunction* clone() const {
		return (ADoubleFunction*)(new SigmaFunction(*this));
	}
};

class SigmaPrimeFunction : public ADoubleFunction {
private:
	double arg0; ///< Sigma's zero
	double kright; ///< right slope
	double kleft; ///< left slope

public:
	virtual ~SigmaPrimeFunction() {};

	SigmaPrimeFunction(double const phi0, double const k1, double const s1)
		: arg0(phi0), kright(k1), kleft(s1) {}

	virtual double operator()(double const arg) const {
		if (arg >= arg0) {
			return kright;
		} else {
			return kleft;
		}
	}

	virtual ADoubleFunction* clone() const {
		return (ADoubleFunction*)(new SigmaPrimeFunction(*this));
	}
};

#endif


