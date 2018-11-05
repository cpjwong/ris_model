import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg
from mpl_toolkits.mplot3d import Axes3D
from random import randint
import warnings
warnings.filterwarnings("ignore")
import sys
def cos(x):
	x = x*np.pi/180
	f = np.cos(x)
	return f
def sin(x):
	x = x*np.pi/180
	f = np.sin(x)
	return f
def lamb1(sigma,omega):
	a = 0.5*(sigma*(1+omega)+1+np.sqrt((-sigma*(1+omega)+1)**2+8*sigma))
	return a
def lamb2(sigma,omega):
	a = 0.5*(sigma*(1+omega)+1-np.sqrt((-sigma*(1+omega)+1)**2+8*sigma))
	return a
def Z(n,lamb1,lamb2):
	f=(lamb1**(n-1)*(lamb2-1))/(lamb2-lamb1)+(lamb2**(n-1)*(1-lamb1))/(lamb2-lamb1)
	return f
for i2 in range(len(sys.argv)):
	if sys.argv[i2]=="-n":
		n=int(sys.argv[i2+1])
	else:
		pass
#Evaluation of partition function Z and probability
sigma = 0.54
omega=0.088
l = 1.54
theta= 180-112
l1 = lamb1(sigma,omega)
l2 = lamb2(sigma,omega)
U = np.zeros((2,2))
U[0][0] = 1
U[1][0] = 1
U[0][1] = 2*sigma
U[1][1] = sigma*(1+omega)
v,w = linalg.eig(U)
w_in = linalg.inv(w)
plt.figure()
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)
plt.gcf().set_size_inches(4,3,forward=True)
plt.subplots_adjust(left=0.18,bottom=0.15)
#n=100
n1 = np.arange(2,n,1)
count = 0
bond_i_a = np.arange(2,n,1)
p_i = []
for i in range(len(bond_i_a)):
	bond_i = bond_i_a[i]
	xl = bond_i-2
	xr = n-1-bond_i
	U_prim = np.zeros((2,2))
	U_prim[0][1] = 2*sigma
	U_prim[1][1] = sigma*(1+omega)
	vd = np.zeros((2,2))
	if xl == 0:
		vd[0][0] = v[0]**xr
		vd[1][1] = v[1]**xr
		Tr = np.mat(w)*np.mat(vd)
		Tr = np.mat(Tr)*np.mat(w_in)
		T= np.mat(U_prim)*np.mat(Tr)
	elif xr == 0:
		vd[0][0] = v[0]**xl
		vd[1][1] = v[1]**xl
		Tl = np.mat(w)*np.mat(vd)
		Tl = np.mat(Tl)*np.mat(w_in)
		T = np.mat(Tl)*np.mat(U_prim)
	else:
		vd[0][0] = v[0]**xr
		vd[1][1] = v[1]**xr
		Tr = np.mat(w)*np.mat(vd)
		Tr = np.mat(Tr)*np.mat(w_in)
		vd[0][0] = v[0]**xl
		vd[1][1] = v[1]**xl
		Tl = np.mat(w)*np.mat(vd)
		Tl = np.mat(Tl)*np.mat(w_in)
		T = np.mat(Tl)*np.mat(U_prim)
		T = np.mat(T)*np.mat(Tr)
	if xl == 0 and xr == 0:
		T = U_prim
	Zt = Z(n,l1,l2)
	row = T[0,:].reshape(1,2).T
	p = float(0.5*(row[0]+row[1]))
	p_i.append(p/Zt)
	#print(p/Zt)
plt.plot(bond_i_a,p_i,'ob-')
st_s = "$n=$"+str(n)
plt.xlabel("$\mathrm{Bond~}i$")
plt.ylabel("$p_{g^+;i}$")
plt.xlim(1,21)
plt.savefig("i_prob2.png",dpi=300)

n1 = np.arange(2,n,1)
count = 0
bond_i_a = np.arange(3,n,1)
p_i_tg = []
for i in range(len(bond_i_a)):
	bond_i = bond_i_a[i]
	xl = bond_i-2
	xr = n-1-bond_i
	U_prim = np.zeros((2,2))
	U_prim[0][1] = 2*sigma
	vd = np.zeros((2,2))
	if xl == 0:
		vd[0][0] = v[0]**xr
		vd[1][1] = v[1]**xr
		Tr = np.mat(w)*np.mat(vd)
		Tr = np.mat(Tr)*np.mat(w_in)
		T= np.mat(U_prim)*np.mat(Tr)
	elif xr == 0:
		vd[0][0] = v[0]**xl
		vd[1][1] = v[1]**xl
		Tl = np.mat(w)*np.mat(vd)
		Tl = np.mat(Tl)*np.mat(w_in)
		T = np.mat(Tl)*np.mat(U_prim)
	else:
		vd[0][0] = v[0]**xr
		vd[1][1] = v[1]**xr
		Tr = np.mat(w)*np.mat(vd)
		Tr = np.mat(Tr)*np.mat(w_in)
		vd[0][0] = v[0]**xl
		vd[1][1] = v[1]**xl
		Tl = np.mat(w)*np.mat(vd)
		Tl = np.mat(Tl)*np.mat(w_in)
		T = np.mat(Tl)*np.mat(U_prim)
		T = np.mat(T)*np.mat(Tr)
	if xl == 0 and xr == 0:
		T = U_prim
	Zt = Z(n,l1,l2)
	row = T[0,:].reshape(1,2).T
	p = float(0.5*(row[0]+row[1]))
	p_i_tg.append(p/Zt)

n1 = np.arange(2,n,1)
count = 0
bond_i_a = np.arange(3,n,1)
p_i_gt = []
for i in range(len(bond_i_a)):
	bond_i = bond_i_a[i]
	xl = bond_i-2
	xr = n-1-bond_i
	U_prim = np.zeros((2,2))
	U_prim[1][0] = 1
	vd = np.zeros((2,2))
	if xl == 0:
		vd[0][0] = v[0]**xr
		vd[1][1] = v[1]**xr
		Tr = np.mat(w)*np.mat(vd)
		Tr = np.mat(Tr)*np.mat(w_in)
		T= np.mat(U_prim)*np.mat(Tr)
	elif xr == 0:
		vd[0][0] = v[0]**xl
		vd[1][1] = v[1]**xl
		Tl = np.mat(w)*np.mat(vd)
		Tl = np.mat(Tl)*np.mat(w_in)
		T = np.mat(Tl)*np.mat(U_prim)
	else:
		vd[0][0] = v[0]**xr
		vd[1][1] = v[1]**xr
		Tr = np.mat(w)*np.mat(vd)
		Tr = np.mat(Tr)*np.mat(w_in)
		vd[0][0] = v[0]**xl
		vd[1][1] = v[1]**xl
		Tl = np.mat(w)*np.mat(vd)
		Tl = np.mat(Tl)*np.mat(w_in)
		T = np.mat(Tl)*np.mat(U_prim)
		T = np.mat(T)*np.mat(Tr)
	if xl == 0 and xr == 0:
		T = U_prim
	Zt = Z(n,l1,l2)
	row = T[0,:].reshape(1,2).T
	p = float(0.5*(row[0]+row[1]))
	p_i_gt.append(p/Zt)

n1 = np.arange(2,n,1)
count = 0
bond_i_a = np.arange(3,n,1)
p_i_gg = []
for i in range(len(bond_i_a)):
	bond_i = bond_i_a[i]
	xl = bond_i-2
	xr = n-1-bond_i
	U_prim = np.zeros((2,2))
	#U_prim[0][0] = 1
	#U_prim[1][0] = 1
	#U_prim[0][1] = 2*sigma
	U_prim[1][1] = sigma*(1+omega)
	vd = np.zeros((2,2))
	if xl == 0:
		vd[0][0] = v[0]**xr
		vd[1][1] = v[1]**xr
		Tr = np.mat(w)*np.mat(vd)
		Tr = np.mat(Tr)*np.mat(w_in)
		T= np.mat(U_prim)*np.mat(Tr)
	elif xr == 0:
		vd[0][0] = v[0]**xl
		vd[1][1] = v[1]**xl
		Tl = np.mat(w)*np.mat(vd)
		Tl = np.mat(Tl)*np.mat(w_in)
		T = np.mat(Tl)*np.mat(U_prim)
	else:
		vd[0][0] = v[0]**xr
		vd[1][1] = v[1]**xr
		Tr = np.mat(w)*np.mat(vd)
		Tr = np.mat(Tr)*np.mat(w_in)
		vd[0][0] = v[0]**xl
		vd[1][1] = v[1]**xl
		Tl = np.mat(w)*np.mat(vd)
		Tl = np.mat(Tl)*np.mat(w_in)
		T = np.mat(Tl)*np.mat(U_prim)
		T = np.mat(T)*np.mat(Tr)
	if xl == 0 and xr == 0:
		T = U_prim
	Zt = Z(n,l1,l2)
	row = T[0,:].reshape(1,2).T
	p = float(0.5*(row[0]+row[1]))
	p_i_gg.append(p/Zt)


n1 = np.arange(2,n,1)
count = 0
bond_i_a = np.arange(3,n,1)
p_i_tt = []
for i in range(len(bond_i_a)):
	bond_i = bond_i_a[i]
	xl = bond_i-2
	xr = n-1-bond_i
	U_prim = np.zeros((2,2))
	U_prim[0][0] = 1
	vd = np.zeros((2,2))
	if xl == 0:
		vd[0][0] = v[0]**xr
		vd[1][1] = v[1]**xr
		Tr = np.mat(w)*np.mat(vd)
		Tr = np.mat(Tr)*np.mat(w_in)
		T= np.mat(U_prim)*np.mat(Tr)
	elif xr == 0:
		vd[0][0] = v[0]**xl
		vd[1][1] = v[1]**xl
		Tl = np.mat(w)*np.mat(vd)
		Tl = np.mat(Tl)*np.mat(w_in)
		T = np.mat(Tl)*np.mat(U_prim)
	else:
		vd[0][0] = v[0]**xr
		vd[1][1] = v[1]**xr
		Tr = np.mat(w)*np.mat(vd)
		Tr = np.mat(Tr)*np.mat(w_in)
		vd[0][0] = v[0]**xl
		vd[1][1] = v[1]**xl
		Tl = np.mat(w)*np.mat(vd)
		Tl = np.mat(Tl)*np.mat(w_in)
		T = np.mat(Tl)*np.mat(U_prim)
		T = np.mat(T)*np.mat(Tr)
	if xl == 0 and xr == 0:
		T = U_prim
	Zt = Z(n,l1,l2)
	row = T[0,:].reshape(1,2).T
	p = float(0.5*(row[0]+row[1]))
	p_i_tt.append(p/Zt)

qs_gt=np.zeros(len(bond_i_a))
qs_gg=np.zeros(len(bond_i_a))
qs_tt=np.zeros(len(bond_i_a))
qs_tg=np.zeros(len(bond_i_a))
for i in range(len(bond_i_a)):
	qs_gt[i]=p_i_gt[i]/(2*p_i[i])
	qs_gg[i]=p_i_gg[i]/(2*p_i[i])
	qs_tt[i]=p_i_tt[i]/(1-2*p_i[i])
	qs_tg[i]=p_i_tg[i]/(1-2*p_i[i])

fig=plt.figure()
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)
plt.gcf().set_size_inches(3,3,forward=True)
plt.gcf().set_size_inches(3,3,forward=True)
plt.subplots_adjust(left=0.12,bottom=0.15)
ax = fig.add_subplot(111, projection='3d')

xf=[]
yf=[]
zf=[]

n1 = np.arange(2,n,1)
lv = [l,0,0]
lo = [0,0,0]

xf.append(lo[0])
yf.append(lo[1])
zf.append(lo[2])

xf.append(lv[0])
yf.append(lv[1])
zf.append(lv[2])

T1 = np.zeros((3,3))
T1[0][0] = cos(theta)
T1[0][1] = sin(theta)
T1[1][0] = sin(theta)
T1[1][1] = -cos(theta)
T1[2][2] = -1
phi = 180
Tt = np.zeros((3,3))
Tt[0][0] = cos(theta)
Tt[0][1] = sin(theta)
Tt[1][0] = -sin(theta)*cos(phi)
Tt[1][1] = cos(theta)*cos(phi)
Tt[1][2] = -sin(phi)
Tt[2][0] = -sin(theta)*sin(phi)
Tt[2][1] = cos(theta)*sin(phi)
Tt[2][2] = cos(phi)
phi = 60
Tg = np.zeros((3,3))
Tg[0][0] = cos(theta)
Tg[0][1] = sin(theta)
Tg[1][0] = -sin(theta)*cos(phi)
Tg[1][1] = cos(theta)*cos(phi)
Tg[1][2] = -sin(phi)
Tg[2][0] = -sin(theta)*sin(phi)
Tg[2][1] = cos(theta)*sin(phi)
Tg[2][2] = cos(phi)
phi = 300
Tgm = np.zeros((3,3))
Tgm[0][0] = cos(theta)
Tgm[0][1] = sin(theta)
Tgm[1][0] = -sin(theta)*cos(phi)
Tgm[1][1] = cos(theta)*cos(phi)
Tgm[1][2] = -sin(phi)
Tgm[2][0] = -sin(theta)*sin(phi)
Tgm[2][1] = cos(theta)*sin(phi)
Tgm[2][2] = cos(phi)
p_sw=np.zeros(len(n1))
for i in range(len(n1)):
	if i == 0 or i == len(n1)-1:
		pn=randint(0,100)
		pg=randint(1,2)
		if pn <= 23 and pg == 1:
			p_sw[i]=1
		elif pn <= 23 and pg == 2:
			p_sw[i]=2
		elif pn > 23:
			p_sw[i]=0
	else:
		pn=randint(0,100)
		pn = pn/100
		pg=randint(1,2)
		if (p_sw[i-1]==1 or p_sw[i-1]==2) and pn>qs_gt[i] and pn<=qs_gt[i]+qs_gg[i] and pg==1:
			p_sw[i]=1
		elif (p_sw[i-1]==1 or p_sw[i-1]==2) and pn>qs_gt[i] and pn<=qs_gt[i]+qs_gg[i] and pg==2:
			p_sw[i]=2
		elif (p_sw[i-1]==1 or p_sw[i-1]==2) and pn<=qs_gt[i]:
			p_sw[i]=0
		elif (p_sw[i-1]==0) and pn>qs_gt[i]+qs_gg[i] and pn<=qs_gt[i]+qs_gg[i]+qs_tt[i]:
			p_sw[i]=0
		elif (p_sw[i-1]==0) and pn>qs_gt[i]+qs_gg[i]+qs_tt[i] and pn<=qs_gt[i]+qs_gg[i]+qs_tt[i]+qs_tg[i] and pg==1:
			p_sw[i]=1
		elif (p_sw[i-1]==0) and pn>qs_gt[i]+qs_gg[i]+qs_tt[i] and pn<=qs_gt[i]+qs_gg[i]+qs_tt[i]+qs_tg[i] and pg==2:
			p_sw[i]=2

l2=np.zeros(3)
l1=np.mat(T1).dot(lv)
l1 = l1[0,:].reshape(1,3).T
l1 = [float(ii) for ii in l1]
xf.append(l1[0]+l)
yf.append(l1[1])
zf.append(l1[2])
l2[0]=l1[0]+l
l2[1]=l1[1]
l2[2]=l1[2]

for i in range(len(p_sw)):
	if i == 0:
		if p_sw[i]==0:
			l1=np.mat(Tt).dot(lv)
		elif p_sw[i]==1:
			l1=np.mat(Tg).dot(lv)
		elif p_sw[i]==2:
			l1=np.mat(Tgm).dot(lv)
		l1 = l1[0,:].reshape(1,3).T
		l1 = [float(ii) for ii in l1]
		l1=np.mat(T1).dot(l1)
		l1 = l1[0,:].reshape(1,3).T
		l1 = [float(ii) for ii in l1]
		xf.append(l1[0]+l2[0])
		yf.append(l1[1]+l2[1])
		zf.append(l1[2]+l2[2])
		l2[0]=l1[0]+l2[0]
		l2[1]=l1[1]+l2[1]
		l2[2]=l1[2]+l2[2]
	else:
		if p_sw[i]==0:
			l1=np.mat(Tt).dot(lv)
		elif p_sw[i]==1:
			l1=np.mat(Tg).dot(lv)
		elif p_sw[i]==2:
			l1=np.mat(Tgm).dot(lv)
		l1 = l1[0,:].reshape(1,3).T
		l1 = [float(ii) for ii in l1]
		for j in range(i):
			if p_sw[i-1-j]==0:
				l1=np.mat(Tt).dot(l1)
			elif p_sw[i-1-j]==1:
				l1=np.mat(Tg).dot(l1)
			elif p_sw[i-1-j]==2:
				l1=np.mat(Tgm).dot(l1)
			l1 = l1[0,:].reshape(1,3).T
			l1 = [float(ii) for ii in l1]
		l1=np.mat(T1).dot(l1)
		l1 = l1[0,:].reshape(1,3).T
		l1 = [float(ii) for ii in l1]
		xf.append(l1[0]+l2[0])
		yf.append(l1[1]+l2[1])
		zf.append(l1[2]+l2[2])
		l2[0]=l1[0]+l2[0]
		l2[1]=l1[1]+l2[1]
		l2[2]=l1[2]+l2[2]			
xcm=sum(xf)/len(xf)
ycm=sum(yf)/len(yf)
zcm=sum(zf)/len(zf)

rgx=0
rgy=0
rgz=0
for i in range(len(xf)):
	rgx=rgx+(xf[i]-xcm)**2
	rgy=rgy+(yf[i]-ycm)**2
	rgz=rgz+(zf[i]-zcm)**2
rgx=np.sqrt(rgx/len(xf))
rgy=np.sqrt(rgy/len(yf))
rgz=np.sqrt(rgz/len(zf))

rg_magn=np.sqrt(rgx**2+rgy**2+rgz**2)
a=-rg_magn
b=rg_magn
plt.xlim(a+xcm,b+xcm)
plt.ylim(a+ycm,b+ycm)
ax.set_zlim(a+zcm,b+zcm)
ax.plot_wireframe(xf,yf,zf)
ax.scatter(xf,yf,zf,'o',c='b')
plt.savefig("3d-plot.png",dpi=300)

filename2="initial2.pdb"
writefile = open(filename2,'w')
for j in range(n):
	if j+1<10:
		writefile.write("ATOM")
		writefile.write("      ")
		writefile.write(str(j+1))
		writefile.write("  ")
		if j==0 or j==n-1:
			writefile.write("C2")
			writefile.write("   ")
		elif j%2==0 and j!=0 and j!=(n-1):
			writefile.write("C1A")
			writefile.write("  ")
		else:
			writefile.write("C1B")
			writefile.write("  ")
		writefile.write("PEE")
		writefile.write("     ")
		writefile.write("1")
		writefile.write("%10s%8s%8s"%(str(round(xf[j],2)),str(round(yf[j],2)),str(round(zf[j],2))))
		writefile.write("  1.00")
		writefile.write("  ")
		writefile.write("0.00")
		writefile.write("\n")
	elif j+1>=10 and j+1<100:
		writefile.write("ATOM")
		writefile.write("     ")
		writefile.write(str(j+1))
		writefile.write("  ")
		if j==0 or j==n-1:
			writefile.write("C2")
			writefile.write("   ")
		elif j%2==0 and j!=0 and j!=(n-1):
			writefile.write("C1A")
			writefile.write("  ")
		else:
			writefile.write("C1B")
			writefile.write("  ")
		writefile.write("PEE")
		writefile.write("     ")
		writefile.write("1")
		writefile.write("%10s%8s%8s"%(str(round(xf[j],2)),str(round(yf[j],2)),str(round(zf[j],2))))
		writefile.write("  1.00")
		writefile.write("  ")
		writefile.write("0.00")
		writefile.write("\n")
	elif j+1>=100:
		writefile.write("ATOM")
		writefile.write("    ")
		writefile.write(str(j+1))
		writefile.write("  ")
		if j==0 or j==n-1:
			writefile.write("C2")
			writefile.write("   ")
		elif j%2==0 and j!=0 and j!=(n-1):
			writefile.write("C1A")
			writefile.write("  ")
		else:
			writefile.write("C1B")
			writefile.write("  ")
		writefile.write("PEE")
		writefile.write("     ")
		writefile.write("1")
		writefile.write("%10s%8s%8s"%(str(round(xf[j],2)),str(round(yf[j],2)),str(round(zf[j],2))))
		writefile.write("  1.00")
		writefile.write("  ")
		writefile.write("0.00")
		writefile.write("\n")
