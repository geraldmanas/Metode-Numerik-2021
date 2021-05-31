#!/usr/bin/env python
# coding: utf-8

# In[15]:


#!/usr/bin/env python
# coding: utf-8

# In[8]:


#===libraries for Calculations===
import sys
import numpy as np
import matplotlib.pyplot as plt
from IPython import get_ipython
from numpy import float32, single, double, array, zeros, diag, diagflat, dot
from pprint import pprint
from scipy.linalg import solve

def modul_2():
    # Metode Setengah Interval
    def setengah_interval(X1,X2):
        print(" =========================================================\n",
              "Akar-Akar Persamaan: Metode Setengah Interval\n",
              "Nama: Gerald Alfa Daud Manas \n",
              "NIM: 26050119140121 \n",
              "Kelas: Oseanografi B \n",
              "=========================================================\n")
        X1 = X1
        X2 = X2
        error = 1
        iterasi = 0
        while(error > 0.001):
            iterasi +=1
            FXi = (float32(-0.1546*X1)**3)+(float32(3.1816*X1)**2)+(float32(-14.097*X1))+(float32(-4.8275))
            FXii = (float32(-0.1546*X2)**3)+(float32(3.1816*X2)**2)+(float32(-14.097*X2))+(float32(-4.8275))
            Xt = (float32(X1 + X2)/2)
            FXt = (float32(-0.1546*Xt)**3)+(float32(3.1816*Xt)**2)+(float32(-14.097*Xt))+(float32(-4.8275))
            if FXi * FXt > 0:
                X1 = Xt
            elif FXi * FXt < 0:
                X2 = Xt
            else:
                print("Akar Penyelesaian: ", Xt)
            if FXt < 0:
                error = FXt * (-1)
            else:
                error = FXt
            if iterasi > 100:
                print("Angka tak hingga")
                break
            print(iterasi, "|", FXi, "|", FXii, "|", Xt, "|", FXt, "|", error)
        print("Jumlah Iterasi: ", iterasi)
        print("Akar persamaan adalah: ", Xt)
        print("Toleransi Error: ", error)
    # Metode Interpolasi Linier
    def interpolasi_linier(X1,X2):
        print(" =========================================================\n",
              "Akar-Akar Persamaan: Metode Interpolasi Linier\n",
              "Nama: Gerald Alfa Daud Manas  \n",
              "NIM: 26050119140121 \n",
              "Kelas: Oseanografi B \n",
              "=========================================================\n")
        X1 = X1
        X2 = X2
        error = 1
        iterasi = 0
        while(error > 0.001):
            iterasi +=1
            FX1 = (float32(-0.1546*X1)**3)+(float32(3.1816*X1)**2)+(float32(-14.097*X1))+(float32(-4.8275))
            FX2 = (float32(-0.1546*X2)**3)+(float32(3.1816*X2)**2)+(float32(-14.097*X2))+(float32(-4.8275))
            Xt = X2 - (single(FX2/(FX2-FX1)))*(X2-X1)
            FXt = (float32(-0.1546*Xt)**3)+(float32(3.1816*Xt)**2)+(float32(-14.097*Xt))+(float32(-4.8275))
            if FXt*FX1 > 0:
                X2 = Xt
                FX2 = FXt
            else:
                X1 = Xt
                FX1 = FXt 
            if FXt < 0:
                error = FXt * (-1)
            else:
                error = FXt
            if iterasi > 100:
                print("Angka tak hingga")
                break
            print(iterasi, "|", FX1, "|", FX2, "|", Xt, "|", FXt, "|", error)
        print("Jumlah Iterasi: ", iterasi)
        print("Akar persamaan adalah: ", Xt)
        print("Toleransi Error: ", error)
    # Metode Metode Newton-Raphson
    def newton_raphson(X1):
        print(" =========================================================\n",
              "Akar-Akar Persamaan: Metode Newton-Raphson\n",
              "Nama: Gerald Alfa Daud Manas \n",
              "NIM: 26050119140121 \n",
              "Kelas: Oseanografi B \n",
              "=========================================================\n")
        X1 = X1
        iterasi = 0
        akar = 1
        while (akar > 0.001):
            iterasi += 1
            Fxn = (float32(-0.1546*X1)**3)+(float32(3.1816*X1)**2)+(float32(-14.097*X1))+(float32(-4.8275))
            Fxxn = (float32(3*X1)**2)+(float32(2*X1))+(float32(1))
            xnp1 = (double(X1 - (Fxn/Fxxn)))
            fxnp1 = (float32(-0.1546*xnp1)**3)+(float32(3.1816*xnp1)**2)+(float32(-14.097*xnp1))+(float32(-4.8275))
            Ea = (single(xnp1-X1)/xnp1)*100
            if Ea < 0.001:
                X1 = xnp1
                akar = Ea*(-1)
            else:
                akar = xnp1
                print("Nilai akar adalah: ", akar)
                print("Nilai error adalah: ", Ea)
            if iterasi > 100:
                break
            print(iterasi, "|", X1, "|", xnp1, "|", akar, "|", Ea)
        print("Jumlah Iterasi: ", iterasi)
        print("Akar persamaan adalah: ", xnp1)
        print("Toleransi Error: ", akar)
    # Metode Secant
    def secant(X1,X2):
        print(" =========================================================\n",
              "Akar-Akar Persamaan: Metode Secant\n",
              "Nama: Gerald Alfa Daud Manas \n",
              "NIM: 26050119140121 \n",
              "Kelas: Oseanografi B \n",
              "=========================================================\n")
        X1 = X1
        X2 = X2
        error = 1
        iterasi = 0
        while(error > 0.001):
            iterasi +=1
            FX1 = (float32(-0.1546*X1)**3)+(float32(3.1816*X1)**2)+(float32(-14.097*X1))+(float32(-4.8275))
            FXmin = (float32(-0.1546*X2)**3)+(float32(3.1816*X2)**2)+(float32(-14.097*X2))+(float32(-4.8275))
            X3 = X1 - (single(FX1))*(single(X1-(X2)))/(single(FX1)-(FXmin))
            FXplus = (float32(-0.1546*X3)**3)+(float32(3.1816*X3)**2)+(float32(-14.097*X3))+(float32(-4.8275))
            if FXplus < 0:
                error = FXplus * (-1)
            else:
                error = FXplus
            if error > 0.001:
                X2 = X1
                X1 = X3
            else:
                print("Selesai")
            if iterasi > 100:
                print("Angka tak hingga")
                break
            print(iterasi, "|", FX1, "|", FXmin, "|", X3, "|", FXplus, "|", error)
        print("Jumlah Iterasi: ", iterasi)
        print("Akar persamaan adalah: ", X3)
        print("Toleransi Error: ", error)
    # Metode Iterasi
    def iterasi(X1):
        print(" =========================================================\n",
              "Akar-Akar Persamaan: Metode Iterasi\n",
              "Nama: Gerald Alfa Daud Manas \n",
              "NIM: 26050119140121 \n",
              "Kelas: Oseanografi B \n",
              "=========================================================\n")
        X1 = X1
        error = 1
        iterasi = 0
        while (error > 0.001):
            iterasi +=1
            Fxn = (float32(-0.1546*X1)**3)+(float32(3.1816*X1)**2)+(float32(-14.097*X1))+(float32(-4.8275))
            X2 = ((float32(-3.1816*X1)**2)+(float32(14.097*X1))+(float32(4.8275)/(float32(-0.1546))))**(float32((0.333334)))
            Ea = (float32((X2-X1)/(X2))*100)
            if Ea < error:
                X1 = X2
                if Ea > 0:
                    error = Ea
                else:
                    error = Ea*(-1)
            else:
                error = Ea
            if iterasi > 100:
                print("Angka tak hingga")
                break
            print(iterasi, "|", X1, "|", X2, "|", Ea, "|", error)
        print("Jumlah Iterasi: ", iterasi)
        print("Akar persamaan adalah: ", X2)
        print("Toleransi Error: ", error)
    print("Kode untuk penyelesaian akar persamaan: \n",
         "1. Metode Setengah Interval \n",
         "2. Metode Interpolasi Linier \n",
         "3. Metode Newton-Raphson \n", 
         "4. Metode Secant \n",
         "5. Metode Iterasi \n")
    setting = int(input("Masukkan kode penyelesaian akar persamaan: "))
    if (setting == 1):
        X1 = (float32(input("Masukkan Nilai Pertama: "))) # Untuk X1 dapat diisi 20
        X2 = (float32(input("Masukkan Nilai Kedua: "))) # Untuk X2 dapat diisi -20
        Cel = setengah_interval(X1,X2)
    elif (setting == 2):
        X1 = (float32(input("Masukkan Nilai Pertama: ")))
        X2 = X1 + 1
        Cel = interpolasi_linier(X1,X2)
    elif (setting == 3):    
        X1 = (float32(input("Masukkan Nilai Pertama: ")))
        Cel = newton_raphson(X1)
    elif (setting == 4):    
        X1 = (float32(input("Masukkan Nilai Pertama: ")))
        X2 = X1 - 1
        Cel = secant(X1,X2)
    elif (setting == 5):    
        X1 = (float32(input("Masukkan Nilai Pertama: ")))
        Cel = iterasi(X1)
    else:
        print("Periksa Kembali Kode yang Diminta!")
        
def modul_3():
    # Metode Gauss
    def Gauss(A, f):
        A = np.array((A), dtype=float)
        f = np.array((f), dtype=float)
        n = len(f)
        for i in range(0, n - 1):  # Looping untuk kolom matriks
            if np.abs(A[i, i]) == 0:
                for k in range(i + 1, n):
                    if np.abs(A[k, i]) > np.abs(A[i, i]):
                        A[[i, k]] = A[[k, i]]  # Tukar antara baris i dan k
                        f[[i, k]] = f[[k, i]]
                        break
            for j in range(i + 1, n):  # Ulangi baris yang ada di bawah diagonal untuk setiap kolom
                m = A[j, i] / A[i, i]
                A[j, :] = A[j, :] - m * A[i, :]
                f[j] = f[j] - m * f[i]
        return A, f
    # Metode Gauss Jordan
    def GaussJordan(a,n):
        #Step1 ===> Looping untuk pengolahan metode Gauss Jordan
        print('==============Mulai Iterasi===============')
        for i in range(n):
            if a[i][i]==0:
                sys.exit('Dibagi dengan angka nol (proses tidak dapat dilanjutkan)')
            for j in range(n):
                if i !=j:
                    ratio=a[j][i]/a[i][i]
                    #print('posisi nol di:[',j,i,']', 'nilai rasio:',ratio)
                    for k in range(n+1):
                        a[j,k]=a[j][k]-ratio*a[i][k]
                    print(a)
                    print(f'============================================')
        # Step 2 ====> Membuat semua variabel(x,y,z,...)==1
        ax=np.zeros((n,n+1))
        for i in range(n):
            for j in range(n+1):
                ax[i,j]=a[i][j]/a[i][i]
        print('===================Akhir Iterasi===================')
        return ax
    # Metode Jacobi
    def Jacobi(A,b, N=10, x=None, info=True):
        print("Hasil perhitungan : ")
        if x is None:
            x = zeros(len(A[0]))
        D = diag(A)
        R = A - diagflat(D)
        for i in range(N):
            x = (b - dot(R,x))/D
            if info:
                pprint(x)
        return x
    # Metode Gauss Seidel
    def GaussSeidel(A,B,x,n):
        L= np.tril(A)
        U=A-L
        for i in range(n):
            x = np.dot(np.linalg.inv(L), B-np.dot(U, x))
            print(x)
        return x
    print("Kode untuk penyelesaian sistem persamaan linier dan matriks: \n",
         "1. Metode Gauss \n",
         "2. Metode Gauss Jordan \n",
         "3. Metode Jacobi \n", 
         "4. Metode Gauss Seidel \n")
    setting = int(input("Masukkan kode penyelesaian sistem persamaan linier dan matriks: "))
    if (setting == 1):
        print(" ================================================ \n",
              "Sistem Persamaan Linier dan Matriks: Metode Gauss \n",
              "Nama: Gerald Alfa Daud Manas \n",
              "NIM: 26050119140121 \n",
              "Kelas: Oseanografi B \n",
              "================================================ \n")
        A = np.array([[2, -6, -7, 4], [6, -7, -7, -1], [-6, 4, 5, -2], [-3, 2, 4, 6]])
        f = np.array([3.5, -0.54, -1, 4.8])
        print('A = \n%s and f = %s \n' % (A, f))
        Cel = Gauss(A, f)
        x = np.linalg.solve(A, f)
        print('Hasil perhitungan dengan metode Gauss adalah x = \n %s \n' % x)
    elif (setting == 2):
        print(" ================================================ \n",
              "Sistem Persamaan Linier dan Matriks: Metode Gauss Jordan \n",
              "Nama: Gerald Alfa Daud Manas \n",
              "NIM: 26050119140121 \n",
              "Kelas: Oseanografi B \n",
              "================================================ \n")
        m = np.array([[2, -6, -7, 4, 3.5],
                      [6, -7, -7, -1, -0.54],
                      [-6, 4, 5, -2, -1],
                      [-3, 2, 4, 6, 4.8]],dtype=float)
        n = 4
        print('Matriks Persamaan') # Menampilkan matrix awal
        print(m)
        m = GaussJordan(m,n) # Menampilkan Hasil
        print(f"""Hasil Pengolahan menggunakan metode Gauss Jordan didapatkan hasil sebagai berikut:
{m} \n""")
    elif (setting == 3):
        print(" ================================================ \n",
              "Sistem Persamaan Linier dan Matriks: Metode Jacobi \n",
              "Nama: Gerald Alfa Daud Manas \n",
              "NIM: 26050119140121 \n",
              "Kelas: Oseanografi B \n",
              "================================================ \n")
        A = np.array([[2, -6, -7, 4], [6, -7, -7, -1], [-6, 4, 5, -2], [-3, 2, 4, 6]])
        b = np.array([3.5, -0.54, -1, 4.8])
        guess = np.array([0,0,0,0])
        Solusi = Jacobi(A,b, N=5, x=guess)
        print("A :")
        pprint(A)
        print("b :")
        pprint(b)
        print("x :")
        pprint(Solusi)
    elif (setting == 4):
        print(" ================================================ \n",
              "Sistem Persamaan Linier dan Matriks: Metode Gauss Seidel \n",
              "Nama: Gerald Alfa Daud Manas \n",
              "NIM: 26050119140121 \n",
              "Kelas: Oseanografi B \n",
              "================================================ \n")
        A = np.array([[2, -6, -7, 4], [6, -7, -7, -1], [-6, 4, 5, -2], [-3, 2, 4, 6]]) # Untuk masuk ke script input, rumus A cara menulisnya : np.array([[2, -6, -7, 4], [6, -7, -7, -1], [-6, 4, 5, -2], [-3, 2, 4, 6]])
        B = np.array([3.5, -0.54, -1, 4.8]) # Untuk masuk ke script input, rumus B cara menulisnya : [hasil nilai D pers.1, hasil nilai D pers.2, hasil nilai D pers.3, hasil nilai D pers.4]
        x = np.array([0,0,0,0]) # Disarankan untuk memasukan nilai 0,0,0,0 untuk nilai awal x
        n = eval(input('Masukan Batas Iterasi:')) # Disarankan batas iterasi dimulai dari 1 sampai 15 agar tidak terlalu banyak
        Cel = GaussSeidel (A,B,x,n)
        print('Hasil Perhitungan Gauss-Seidel: \n', solve(A,B))
    else:
        print("Periksa Kembali Kode yang Diminta!")
        
def modul_4():
    #Trapesium 1 Pias
    def trapesium_1pias():
        print(" ================================================ \n",
              "Integrasi Numerik: Metode Trapesium 1 Pias \n",
              "Nama: Gerald Alfa Daud Manas \n",
              "NIM: 26050119140121 \n",
              "Kelas: Oseanografi B \n",
              "================================================ \n")
        x = np.linspace(-10,10,100)
        y = 3*(x**3) + 3*(x**2) + 3
        plt.plot(x,y)
        x1 = 0
        x2 = 1
        fx1 = 3*x1**3
        fx2 = 3*x2**2
        plt.fill_between([x1,x2], [fx1,fx2])
        plt.xlim([-1.5,1.5]); plt.ylim([0,4]);
        plt.title('Trapesium 1 pias')
        plt.savefig('C:/Users/Hp/Documents/Tugas Metnum/image/Trapesium 1 Pias.png')
        L = 0.5*(fx2 + fx1)*(x2 - x1)
        print("Luas dengan metode trapesium 1 pias:", L)
    #Trapesium Banyak Pias
    def trapesium_banyakpias(f,a,b,N):
        x = np.linspace(a,b,N+1)
        y = f(x)
        y_kanan = y[1:]
        y_kiri = y[:-1]
        dx = (b - a)/N
        T = (dx/2) * np.sum(y_kanan + y_kiri)
        return T
    f = lambda x : 3*(x**3) + 3*(x**2) + 3
    a = 0
    b = 10
    N = 3
    x = np.linspace(a,b,N+1)
    y = f(x)
    X = np.linspace(a,b+1,N)
    Y = f(X)
    #Metode Simpson 1/3
    def simpson1per3(x0,xn,n):
        f = lambda x: 3*(x**3) + 3*(x**2) + 3
        h = (xn - x0) / n
        integral = f(x0) + f(xn)
        for i in range(1,n):
            k = x0 + i*h
            if i%2 == 0:
                integral = integral + 2 * f(k)
            else:
                integral = integral + 4 * f(k)
        integral = integral * h / 3
        return integral
    #Simpson 3/8
    def simpson3per8(x0,xn,n):
        print(" ================================================ \n",
              "Integrasi Numerik: Metode Simpson 3/8 \n",
              "Nama: Gerald Alfa Daud Manas \n",
              "NIM: 26050119140121 \n",
              "Kelas: Oseanografi B \n",
              "================================================ \n")
        f = lambda x: 3*(x**3)+3*(x**2)
        h = (xn - x0) / n
        integral = f(x0) + f(xn)
        for i in range(1,n):
            k = x0 + i*h
            if i%2 == 0:
                integral = integral + 3 * f(k)
            else:
                integral = integral + 3 * f(k)
        integral = integral * 3 * h / 8
        return integral
    print("Kode untuk penyelesaian integrasi numerik: \n",
         "1. Metode Trapesium 1 Pias \n",
         "2. Metode Trapesium Banyak Pias \n",
         "3. Metode Simpson 1/3 \n",
         "4. Metode Simpson 3/8 \n")
    setting = int(input("Masukkan kode penyelesaian integrasi numerik: "))
    if (setting == 1):
        Cel = trapesium_1pias()
    elif (setting == 2):
        print(" ================================================ \n",
             "Integrasi Numerik: Metode Trapesium Banyak Pias \n",
             "Nama  : Gerald Alfa Daud Manas \n",
             "NIM   : 26050119140121 \n",
             "Kelas : Oseanografi B \n",
             "================================================ \n")
        plt.plot(X,Y)
        for i in range(N):
            xs = [x[i],x[i],x[i+1],x[i+1]]
            ys = [0,f(x[i]),f(x[i+1]),0]
            plt.fill(xs,ys, 'b', edgecolor='b',alpha=0.2)
        plt.title('Trapesium banyak pias, N = {}'.format(N))
        plt.savefig('C:/Users/Hp/Documents/Tugas Metnum/image/Trapesium Banyak Pias.png')
        L = trapesium_banyakpias(f,a,b,N)
        print(L)
    elif (setting == 3):
        print(" ================================================ \n",
             "Integrasi Numerik: Metode Simpson 1/3 \n",
             "Nama  : Gerald Alfa Daud Manas \n",
             "NIM   : 26050119140121 \n",
             "Kelas : Oseanografi B \n",
             "================================================ \n")
        x1 = float(input("Batas bawah (a): "))
        x2 = float(input("Batas bawah (b): "))
        hasil = simpson1per3(x1, x2, 3)
        print("Nilai integral metode Simpson 1/3:", hasil)
    elif (setting == 4):
        x1 = float(input("Batas bawah (x1): "))
        x2 = float(input("Batas atas (x2): "))  #Isi dengan nilai 10 karna sesuai dalam tugas h+6
        hasil = simpson3per8(x1, x2, 3)
        print("Nilai integral metode Simpson 3/8:", hasil)
    else:
        print("Periksa Kembali Kode yang Diminta!")
    
def modul_5():
    # Metode Euler
    def euler():
        Nama = (input("Masukkan Nama: "))
        Nim = (input("Masukkan NIM: "))
        Kelas = (input("Masukkan Kelas: "))
        plt.style.use('seaborn-poster')
        ipy = get_ipython()
        if ipy is not None:
            ipy.run_line_magic('matplotlib', 'inline')
        h = float(input("Masukkan nilai h: ")) #Banyak langkah yang akan digunakan (Ose B = 0.05)
        x0 = float(input("Masukkan nilai x awal: ")) #Masukkan kondisi awal 0
        xn = float(input("Masukkan nilai x akhir: ")) #Untuk x akhir ketentuan 2 nim akhir dibalik (nilainya 3,0)
        x = np.arange(x0, xn + h, h) #Grid
        y0 = float(input("Masukkan nilai y awal: ")) #Masukkan kondisi awal 0
        G = 2*(x**3) + 9*(x**2) + 3*(x) + 3 #Fungsi
        f = lambda x, y: 2*(x**3) + 9*(x**2) + 3*(x) + 3 #Persamaan Diferensial Biasa
        y = np.zeros(len(x))
        y[0] = y0
        for i in range(0, len(x) - 1):
            y[i + 1] = y[i] + h*f(x[i], y[i])
        print(y)
        error = G-y
        print(error)
        judul = ("\n Grafik Pendekatan Persamaan Differensial Biasa Dengan Metode Euler \n\n %s_%s_%s \n" % (Nama,Nim,Kelas))
        plt.figure(figsize = (12, 10))
        plt.plot(x, y, 'b--', color='r', label='Hasil Pendekatan')
        plt.plot(x, G, '-g', color='y', label='Hasil Analitik')
        plt.title(judul) #Judul Plot
        plt.xlabel('x') #Label u/ sumbu x
        plt.ylabel('F(x)') #Label u/ sumbu y
        plt.grid() #u/ menampilkan grid
        plt.legend(loc='lower right')
        plt.savefig('C:/Users/Hp/Documents/Tugas Metnum/image/Euler.png')
    # Metode Heun
    def heun():
        Nama = (input("Masukkan Nama: "))
        Nim = (input("Masukkan NIM: "))
        Kelas = (input("Masukkan Kelas: "))
        plt.style.use('seaborn-poster')
        ipy = get_ipython()
        if ipy is not None:
            ipy.run_line_magic('matplotlib', 'inline')
        h = float(input("Masukkan nilai h: ")) #Banyak langkah yang akan digunakan (Ose B = 0.05)
        x0 = float(input("Masukkan nilai x awal: ")) #Masukkan kondisi awal 0
        xn = float(input("Masukkan nilai x akhir: ")) #Untuk x akhir ketentuan 2 nim akhir dibalik (nilainya 3,0)
        x = np.arange(x0, xn + h, h) #Grid
        y0 = float(input("Masukkan nilai y awal: ")) #Masukkan kondisi awal 0
        G = 2*(x**3) + 9*(x**2) + 3*(x) + 3 #Fungsi
        f = lambda x, y: 2*(x**3) + 9*(x**2) + 3*(x) + 3 #Persamaan Diferensial Biasa
        y = np.zeros(len(x))
        y[0] = y0

        for i in range(0, len(x) - 1):
            k1 = f(x[i], y[i])
            k2 = f((x[i] + h), (y[i] + (h*k1)))
            y[i+1] = y[i] + (0.5 * h * (k1 + k2))
        print(y)
        error = G-y
        print(error)
        judul = ("\n Grafik Pendekatan Persamaan Differensial Biasa Dengan Metode Heun \n\n %s_%s_%s \n" % (Nama,Nim,Kelas))
        plt.figure(figsize = (12, 10))
        plt.plot(x, y, 'b--', color='r', label='Hasil Pendekatan')
        plt.plot(x, G, '-g', color='y', label='Hasil Analitik')
        plt.title(judul) #Judul Plot
        plt.xlabel('x') #Label u/ sumbu x
        plt.ylabel('F(x)') #Label u/ sumbu y
        plt.grid() #u/ menampilkan grid
        plt.legend(loc='lower right')
        plt.savefig('C:/Users/Hp/Documents/Tugas Metnum/image/Heun.png')
    print("Kode untuk penyelesaian persamaan differensial biasa: \n",
         "1. Metode Euler \n",
         "2. Metode Heun \n")
    setting = int(input("Masukkan kode penyelesaian persamaan differensial biasa: "))
    if (setting == 1):
        Cel = euler()
    elif (setting == 2):
        Cel = heun()
    else:
        print("Lord Said: revisi kode yang diminta")
print("Kode untuk materi metode numerik: \n",
      "1. Modul 2 Akar-akar Persamaan \n",
      "2. Modul 3 Sistem Persamaan Linier dan Matriks \n",
      "3. Modul 4 Integrasi Numerik \n",
      "4. Modul 5 Persamaan Persamaan Differensial Biasa \n")
setting = int(input("Masukkan kode untuk materi metode numerik: "))
if (setting == 1):
    Cel = modul_2()
elif (setting == 2):
    Cel = modul_3()
elif (setting == 3):
    Cel = modul_4()
elif (setting == 4):
    Cel = modul_5()
else:
    print("Lord Said: Kalo tahu salah coba intropeksi diri, salahnya dimana, kodenya jangan template!")


# In[ ]:




