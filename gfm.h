/* Tumour cord growth model */

/*
 * Copyright (C) 2007 Sergey Astanin, Luigi Preziosi, Andrea Tosin
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

#ifndef GFM_H
#define GFM_H

#include "function.h"
#include "amesh2d.h"
#include "newton.h"

#include <sstream>
#include <memory>
using std::auto_ptr;
using std::ostringstream;

class ASigmaFactory {
public:
	virtual ~ASigmaFactory() {}
	virtual ADoubleFunction* build_sigma() = 0;
	virtual ADoubleFunction* build_sigma_prime() = 0;
};

/// piece linear tumour sigma
class TumourSigma : public ASigmaFactory {
	double phi0;
	double kright;
	double kleft;
public:
	virtual ~TumourSigma() {}
	// may throw MeshException
	TumourSigma(AMesh2D const& m) : phi0(m.get_attr("phi0")),
		kright(m.get_attr("tk1")), kleft(m.get_attr("ts1")) {
	}
	virtual ADoubleFunction* build_sigma() {
		return new SigmaFunction(phi0,kright,kleft);
	}
	virtual ADoubleFunction* build_sigma_prime() {
		return new SigmaPrimeFunction(phi0,kright,kleft);
	}
};

/// piece linear host sigma
class HostSigma : public ASigmaFactory {
	double phi0;
	double kright;
	double kleft;
public:
	virtual ~HostSigma() {}
	// may throw MeshException
	HostSigma(AMesh2D const& m) : phi0(m.get_attr("phi0")),
		kright(m.get_attr("hk1")), kleft(m.get_attr("hs1")) {
	}
	virtual ADoubleFunction* build_sigma() {
		return new SigmaFunction(phi0,kright,kleft);
	}
	virtual ADoubleFunction* build_sigma_prime() {
		return new SigmaPrimeFunction(phi0,kright,kleft);
	}
};

/// converts fluid A to fluid B: phi_A*Sigma(phi_A)=phi_B*Sigma(phi_B)
class GFM_A2B_Converter : public ADoubleFunction {
private:
	// TODO: make them const auto_ptr<>
	const auto_ptr<ADoubleFunction> s_a;
	const auto_ptr<ADoubleFunction> s_b;
	const auto_ptr<ADoubleFunction> sp_b;
	// auxilary functions
	LinearFunction x;
	ProductFunction xs_b; // x*Sigma_b(x)
	ProductFunction xsp_b; // x*Sigma_prime_b(x)
	SumFunction Fprime;
public:
	virtual ~GFM_A2B_Converter() {}

	// build f:  F=x*sigma_b(x)-phi_t*sigma_a(phi_a)
	//  and f':  Fprime=sigma_b(x)+x*sigma_prime_b(x)
	// and use them to find x(phi_a) in operator(phi_a) (Newton method)
	GFM_A2B_Converter(ADoubleFunction * sigma_a,
			ADoubleFunction * sigma_b,
			ADoubleFunction * sigma_prime_b) :
		s_a(sigma_a),
		s_b(sigma_b),
		sp_b(sigma_prime_b),
       		x(1,0), xs_b(x,*s_b), xsp_b(x,*sp_b), Fprime(*s_b,xsp_b) {}

	virtual double
	operator()(double const arg) const {
		ConstFunction T( -arg * s_a->eval(arg) ); // -arg*Sigma_a(arg)
		SumFunction F(xs_b,T); // x*Sigma_b(x)-arg*Sigma_a(arg)
		NewtonRootFinder rf(F,Fprime,arg,1e-10);
		int max_iter=1000;
		int i=0;
		while ((i<max_iter) && (!rf.is_found())) {
			rf.iterate();
			++i;
		}
		if (rf.is_found()) {
			return rf.get_root();
		} else {
			ostringstream ss;
			ss << "GFM_A2B_Converter::operator(): " << "A=" << arg
				<< ", B was not found in" << i << " iterations";
			throw FunctionException(ss.str());
		}
	}
	virtual ADoubleFunction*
	clone() const {
		return new GFM_A2B_Converter(&*s_a,&*s_b,&*sp_b);
	}
};

/// Ghost Fluid Method: construct ghost fluids and reconstruct real fluids
class AGFMSplitter {
public:
	virtual ~AGFMSplitter() {}
	/// split mesh variable fid_real into fid_pos and fid_neg
	virtual void
	split(AMesh2D& m, string const fid_real, string const fid_phase,
		string const fid_pos, string const fid_neg) const = 0;
	/// reconstruct mesh variable fid_real from fid_pos and fid_neg
	virtual void
	merge(AMesh2D& m, string const fid_real, string const fid_phase,
		string const fid_pos, string const fid_neg) const = 0;
};

class TumourHostSplitter : public AGFMSplitter {
private:
	TumourSigma tumour;
	HostSigma host;
	GFM_A2B_Converter t2h;
	GFM_A2B_Converter h2t;
public:
	virtual ~TumourHostSplitter() {}
	TumourHostSplitter(AMesh2D const& m) :
		tumour(m), host(m),
		t2h(tumour.build_sigma(),
			host.build_sigma(),host.build_sigma_prime()),
		h2t(host.build_sigma(),
			tumour.build_sigma(),tumour.build_sigma_prime()) {}
	/// split mesh variable fid_real into fid_pos and fid_neg
	virtual void
	split(AMesh2D& m, string const fid_real, string const fid_phase,
		string const fid_pos, string const fid_neg) const {
		try {
		m.remove_function_ifdef(fid_pos);
		m.remove_function_ifdef(fid_neg);
		m.add_function(fid_pos);
		m.add_function(fid_neg);
		for (int i=0; i<m.get_xdim(); ++i) {
			for (int j=0; j<m.get_ydim(); ++j) {
				double psi=m.get(fid_phase,i,j);
				double phi=m.get(fid_real,i,j);
				double phi_pos, phi_neg;
				if (psi > 0) { // `positive' fluid region
					phi_pos=phi;
					phi_neg=t2h(phi);
				} else { // `negative' fluid region
					phi_pos=h2t(phi);
					phi_neg=phi;
				}
				m.set(fid_pos,i,j,phi_pos);
				m.set(fid_neg,i,j,phi_neg);
			}
		}
		} catch (MeshException& e) {
			ostringstream ss;
			ss << "TumourHostSplitter::split: " << e.what();
			throw MeshException(ss.str());
		}
	}
	/// reconstruct mesh variable fid_real from fid_pos and fid_neg
	virtual void
	merge(AMesh2D& m, string const fid_real, string const fid_phase,
		string const fid_pos, string const fid_neg) const {
		try {
		m.remove_function_ifdef(fid_real);
		m.add_function(fid_real);
		for (int i=0; i<m.get_xdim(); ++i) {
			for (int j=0; j<m.get_ydim(); ++j) {
				double psi=m.get(fid_phase,i,j);
				double phi_pos=m.get(fid_pos,i,j);
				double phi_neg=m.get(fid_neg,i,j);
				double phi_real;
				if (psi > 0) { // `positive' fluid region
					phi_real=phi_pos;
				} else { // `negative' fluid region
					phi_real=phi_neg;
				}
				m.set(fid_real,i,j,phi_real);
			}
		}
		} catch (MeshException& e) {
			ostringstream ss;
			ss << "TumourHostSplitter::merge: " << e.what();
			throw MeshException(ss.str());
		}
	}
};

#endif

