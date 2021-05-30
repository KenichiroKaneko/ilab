import matplotlib.pyplot as plt


r= []
z= []

with open("python/SENPOS1.txt", mode="r", encoding="utf-8") as f:
    for row in f:
        s = row.strip().split(" ")
        if len(s) > 1:
            r.append(float(s[0][:-4])/10)
            if s[1] != "0":
                z.append(float(s[1][:-4])/10)   
            else:
                z.append(0)


print(r, z)
print(len(r), len(z))


zmin= -9.985000000000001e-01
zmax= 9.985000000000001e-01
rmin= 1.081500000000000e-01
rmax= 6.940000000000000e-01
delr= 9.747920133111480e-04
delz= 9.827755905511811e-04

matrix_r = []
matrix_z = []

for rr,zz in zip(r,z):
    print(rr,zz)
    rr -= rmin
    if rr%delr < (delr/2):
        matrix_r.append(int(rr//delr))
    else:
        matrix_r.append(int(rr//delr + 1))
    
    if zz%delz < (delz/2):
        matrix_z.append(int(zz//delz + 1017))
    else:
        matrix_z.append(int(zz//delz + 1018))
    
with open('senpos1rz.txt', mode="w", encoding="utf-8")as f:
    for rr, zz in zip(matrix_r,matrix_z):
        f.writelines("{0} {1}\n".format(rr,zz))

# print(matrix_r,matrix_z)
# plt.figure()
# plt.plot(matrix_r, matrix_z, "o")
# plt.show()