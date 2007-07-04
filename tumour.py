#!/usr/bin/env python

xsize = 1.0 # domain length
ysize = 1.0 # domain width
nx = 80 # elements along x
ny = 80 # elements along y
dx = xsize/nx
dy = ysize/ny

r0 = 0.2 # initial radius of the tumour

t=700.0 # simulation time
dt=0.5 # time step
dump_dt=100.0 # dump period

phi0=0.75 # stress-free density
mu=0.01 # tissue relaxation rate
theta=0.15 # base (survival) consumption
epsilon=0.8 # death rate

c_in=1.0 # oxygen in vessel
alpha=200.0 # rate of oxygen consumption

from fipy.meshes.grid2D import Grid2D
from fipy.variables.cellVariable import CellVariable
from fipy.boundaryConditions.fixedValue import FixedValue
from fipy.boundaryConditions.fixedFlux import FixedFlux
from fipy.terms.implicitDiffusionTerm import ImplicitDiffusionTerm
from fipy.terms.implicitSourceTerm import ImplicitSourceTerm
from fipy.terms.transientTerm import TransientTerm
from fipy.tools import numerix
from fipy.models.levelSet.distanceFunction.distanceVariable \
						import DistanceVariable
from fipy.models.levelSet.advection.advectionEquation \
						import buildAdvectionEquation

def pos(x):
	if x > 0:
		return x
	else:
		return 0.0

def sigma(phi):
	return phi-phi0

def sigma_prime(phi):
	return 1.0

def simulation():
	mesh=Grid2D(nx=nx,ny=ny,dx=dx,dy=dy)
	phi=CellVariable(mesh=mesh,name='cell packing density',value=phi0)
	c=CellVariable(mesh=mesh,name='oxygen',value=c_in)
	psi=DistanceVariable(mesh=mesh,name='level set variable',value=1.0)
	psi0=r0-numerix.sqrt((mesh.getCellCenters()[:,0])**2 + \
			(mesh.getCellCenters()[:,1])**2)
	psi.setValue(psi0)
	# psi.calcDistanceFunction()
	Hpsi=CellVariable(mesh=mesh,name='tumour mask H(psi)',value=0.0)

	# boundary conditions, zero flux is assumed if not specified
	phiBCs=(FixedValue(faces=mesh.getFacesRight(),value=phi0),
		FixedValue(faces=mesh.getFacesTop(),value=phi0))
	cBCs = (FixedValue(faces=mesh.getFacesLeft(),value=c_in))

	# model equations
	phiEq = TransientTerm() == ImplicitDiffusionTerm(coeff=
				mu*(phi*sigma(phi)+phi*phi*sigma_prime(phi)) ) \
			+ ImplicitSourceTerm( phi*pos((1-phi)*c-theta)*Hpsi
				- epsilon*phi*pos(theta-(1-phi)*c)*Hpsi )
	cEq = TransientTerm() == ImplicitDiffusionTerm( coeff=1.0 ) \
			+ ImplicitSourceTerm( -alpha*phi*(1-phi)*c*Hpsi )
	gradpsi=psi.getGrad()
	psiEq = buildAdvectionEquation(advectionCoeff=
			-mu*(phi*sigma(phi)+phi*phi*sigma_prime(phi))
			*phi.getGrad().dot(gradpsi)/
			numerix.sqrtDot(gradpsi,gradpsi))

	# solve the equations
	steps=int(t/dt)
	lastdump=0.0
	for s in range(steps):
		ct=s*dt;
		print ct,' max(phi)=%0.4f'%max(phi),' min(phi)=%0.4f'%min(phi)
		print ct,' max(c)  =%0.4f'%max(c),' min(c)  =%0.4f'%min(c)
		print ct,' max(psi)=%0.4f'%max(psi),' min(psi)=%0.4f'%min(psi)

		# update H(psi)
		Hpsi.setValue(1.0,where=(psi>0))
		Hpsi.setValue(0.0,where=(psi<0))

		phiEq.solve(var=phi,boundaryConditions=phiBCs,dt=dt)

		cEq.solve(var=c,boundaryConditions=cBCs,dt=dt)

		psiEq.solve(psi,dt=dt)

		if (ct-lastdump) >= (dump_dt-1e-3*dt):
			dump2gp(mesh,ct,phi,c,psi,Hpsi)
			lastdump=ct

# saves the variables in a format suitable for cord plotting scripts
def dump2gp(mesh, t, phi, c, psi, Hpsi):
	fname="dmp%09d.gp"%(int(t*10))
	f=open(fname,'w')
	f.write("#time: %e\n"%t)
	f.write("#dimensions: %d %d\n"%(ny,nx))
	f.write("#x: ")
	ccs=mesh.getCellCenters()
	xs = list(set(map(lambda cc: cc[0], ccs)))
	ys = list(set(map(lambda cc: cc[1], ccs)))
	for xp in sorted(xs):
		f.write("%f "%xp)
	f.write("\n")
	f.write("#y: ")
	for xp in sorted(ys):
		f.write("%f "%xp)
	f.write("\n")
	f.write("#x #y #c #phi #psi #phi_growth\n")
	phis=phi.getValue()
	cs=c.getValue()
	psis=psi.getValue()
	hs=Hpsi.getValue()
	points=map(lambda e:(e[0][0],e[0][1],e[1],e[2],e[3],e[4]), \
						zip(ccs,cs,phis,psis,hs))
	last=points[0][1]
	for p in points:
		if p[1] > last:
			f.write("\n")
			last=p[1]
		f.write(reduce(lambda x, y: str(x)+' '+str(y),p))
		f.write("\n")
	f.close

if __name__ == '__main__':
	simulation()

