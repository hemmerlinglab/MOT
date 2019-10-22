import numpy as np
import matplotlib.pyplot as plt
import magpylib as magpy

# magpylib standard units are mm and mT

print('building sources...',end='')
rmax = 14.605 # mm
rmin = 6.985 # mm
rstp = (rmin+rmax)/5 # mm

s1u1 = magpy.source.current.Circular(curr=6,dim=2*(rmin))
s2u1 = magpy.source.current.Circular(curr=6,dim=2*(rmin+rstp))
s3u1 = magpy.source.current.Circular(curr=6,dim=2*(rmin+2*rstp))
s4u1 = magpy.source.current.Circular(curr=6,dim=2*(rmin+3*rstp))
s5u1 = magpy.source.current.Circular(curr=6,dim=2*(rmin+4*rstp))
s6u1 = magpy.source.current.Circular(curr=6,dim=2*(rmin+5*rstp))

s1u2 = magpy.source.current.Circular(curr=6,dim=2*(rmin))
s2u2 = magpy.source.current.Circular(curr=6,dim=2*(rmin+rstp))
s3u2 = magpy.source.current.Circular(curr=6,dim=2*(rmin+2*rstp))
s4u2 = magpy.source.current.Circular(curr=6,dim=2*(rmin+3*rstp))
s5u2 = magpy.source.current.Circular(curr=6,dim=2*(rmin+4*rstp))
s6u2 = magpy.source.current.Circular(curr=6,dim=2*(rmin+5*rstp))

s1d1 = magpy.source.current.Circular(curr=-6,dim=2*(rmin))
s2d1 = magpy.source.current.Circular(curr=-6,dim=2*(rmin+rstp))
s3d1 = magpy.source.current.Circular(curr=-6,dim=2*(rmin+2*rstp))
s4d1 = magpy.source.current.Circular(curr=-6,dim=2*(rmin+3*rstp))
s5d1 = magpy.source.current.Circular(curr=-6,dim=2*(rmin+4*rstp))
s6d1 = magpy.source.current.Circular(curr=-6,dim=2*(rmin+5*rstp))

s1d2 = magpy.source.current.Circular(curr=-6,dim=2*(rmin))
s2d2 = magpy.source.current.Circular(curr=-6,dim=2*(rmin+rstp))
s3d2 = magpy.source.current.Circular(curr=-6,dim=2*(rmin+2*rstp))
s4d2 = magpy.source.current.Circular(curr=-6,dim=2*(rmin+3*rstp))
s5d2 = magpy.source.current.Circular(curr=-6,dim=2*(rmin+4*rstp))
s6d2 = magpy.source.current.Circular(curr=-6,dim=2*(rmin+5*rstp))
print('complete')

print('moving sources...',end='')
d1 = 6.55 # mm
d2 = .91 # mm

s1u1.move([0,0,d1+d2])
s2u1.move([0,0,d1+d2])
s3u1.move([0,0,d1+d2])
s4u1.move([0,0,d1+d2])
s5u1.move([0,0,d1+d2])
s6u1.move([0,0,d1+d2])

s1u2.move([0,0,d1-d2])
s2u2.move([0,0,d1-d2])
s3u2.move([0,0,d1-d2])
s4u2.move([0,0,d1-d2])
s5u2.move([0,0,d1-d2])
s6u2.move([0,0,d1-d2])

s1d1.move([0,0,-d1-d2])
s2d1.move([0,0,-d1-d2])
s3d1.move([0,0,-d1-d2])
s4d1.move([0,0,-d1-d2])
s5d1.move([0,0,-d1-d2])
s6d1.move([0,0,-d1-d2])

s1d2.move([0,0,-d1+d2])
s2d2.move([0,0,-d1+d2])
s3d2.move([0,0,-d1+d2])
s4d2.move([0,0,-d1+d2])
s5d2.move([0,0,-d1+d2])
s6d2.move([0,0,-d1+d2])
print('complete')

print('building collection...',end='')
Cmot = magpy.Collection(s1u1,s2u1,s3u1,s4u1,s5u1,s6u1,
					 s1u2,s2u2,s3u2,s4u2,s5u2,s6u2,
					 s1d1,s2d1,s3d1,s4d1,s5d1,s6d1,
					 s1d2,s2d2,s3d2,s4d2,s5d2,s6d2)

Cmot.displaySystem()
print('complete')

xmin = -2*rmax
xmax = 2*rmax
nx = 100
zmin = -2*d1
zmax = 2*d1
nz = 100

xs = np.linspace(xmin,xmax,nx)
zs = np.linspace(zmin,zmax,nz)

Z = [[0,0,z] for z in zs]
X = [[x,0,0] for x in xs]

print('getting z sweep...',end='')
plt.figure()
plt.subplot(211)
plt.title('Z Magnetic Field and dB/dz on Z Axis')
BZ = Cmot.getBsweep(Z) # mT
plt.plot(zs/10,BZ[:,2]*10)
plt.xlabel('z (cm)')
plt.ylabel('Bz (Gs)')
plt.subplot(212)
GZ = np.multiply(np.gradient(BZ[:,2]*10),(10*nz/(zmax-zmin))) # Gs/cm
plt.plot(zs/10,GZ)
plt.xlabel('z (cm)')
plt.ylabel('dBz/dz (Gs/cm)')
print('complete')

print('getting x sweep...',end='')
plt.figure()
plt.title('X Magnetic Field on X Axis')
BX = Cmot.getBsweep(X) # mT
plt.plot(xs/10,BX[:,0]*10)
plt.xlabel('x (cm)')
plt.ylabel('Bx (Gs)')
print('complete')

print('getting xz slice...',end='')
plt.figure()
plt.title('XZ Magnetic Field in XZ Plane')
Bslice = np.array([[Cmot.getB([x,0,z]) for x in xs] for z in zs])
#print(Bslice)
Xm,Zm = np.meshgrid(xs,zs)
U,V = Bslice[:,:,0], Bslice[:,:,2]
plt.pcolor(xs,zs,np.log(U**2+V**2))
plt.streamplot(Xm,Zm,U,V,color='w',density=2)
plt.margins(0,0)
plt.xlabel('x (mm)')
plt.ylabel('z (mm)')
print('complete')

plt.show()