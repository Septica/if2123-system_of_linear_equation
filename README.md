# Tugas Besar 1 IF2123 Aljabar Geometri
## Aplikasi Aljabar Lanjar pada Metode Numerik
### Semester I Tahun 2017/2018

#### DESKRIPSI UMUM TUGAS BESAR
Tugas besar ini adalah membuat program menghitung solusi Sistem Persamaan Linier (SPL) secara numerik dalam bahasa pemrograman Java dengan menggunakan metode eliminasi Gauss dan/atau Gauss-Jordan. SPL dapat memiliki solusi unik, banyak solusi, atau solusi tidak ada. SPL juga digunakan dalam menentukan persamaan polinom interpolasi.

Karena perhitungan menggunakan representasi bilangan titik-kambang (floating point) di dalam komputer, maka untuk meminumkan galat perhitungan, digunakan strategi pivoting dalam memilih baris yang dijadikan basis dalam operasi baris elementer. Bahasa Java digunakan sebagai bahan belajar penggunaan bahasa pemrograman selain C dan Pascal yang sudah digunakan selama ini.

#### SPESIFIKASI UMUM
1. Program harus dapat menerima input data dari
  * Papan ketik
  * File
2. Keluaran program harus dapat ditampilkan ke:
  * Layar monitor
  * Simpan ke dalam arsip
  
#### SPESIFIKASI MATERI
SPL dapat diselesaikan secara numerik dengan metode eliminasi Gauss dan metode eliminasi Gauss-Jordan. Di dalam kedua metode tersebut diterapkan tatancang pemorosan (pivoting) untuk mengurangi galat pembulatan.

Program harus dapat menangani kasus-kasus sebagai berikut:  
a) SPL memiliki solusi unik, tampilkan solusinya  
b) SPL memiliki solusi tak terbatas, tampilkan solusinya dalam bentuk parameter  
c) SPL tidak memiliki solusi, tuliskan tidak ada solusinya.  
