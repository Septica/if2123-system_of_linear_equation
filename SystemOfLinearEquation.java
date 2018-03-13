public class SystemOfLinearEquation {
    private static final double EPSILON = 1e-8;

    private final int M;
    private final int N;

    private double[][] a;

    private boolean echelon;
    private boolean reduced;

    // Constructor of system of linear equation
    public SystemOfLinearEquation(double[][] A, double[] b) {
        // Set matrix size
        M = A.length;
        N = b.length;

        // Build augmented matrix
        a = new double[M][N+1];
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {
                a[i][j] = A[i][j];
            }
        }

        for (int i = 0; i < M; i++) {
            a[i][N] = b[i];
        }
        // Set initial state
        echelon = false;
        reduced = false;
    }

    public SystemOfLinearEquation(double[] xs, double[] fs) {
        // Set matrix size
        M = xs.length;
        N = fs.length;

        // Build augmented matrix
        a = new double[M][N+1];
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {
                a[i][j] = Math.pow(xs[i], j);
            }
        }

        for (int i = 0; i < M; i++) {
            a[i][N] = fs[i];
        }

        // Set initial state
        echelon = false;
        reduced = false;
    }

    // Interchanging row1 and row2 of a matrix
    private void switchRow(int row1, int row2) {
        double[] temp = a[row1];
        a[row1] = a[row2];
        a[row2] = temp;
    }

    // Adding k multiple of row2 to row1
    private void addRow(int row1, int row2, double k) {
        for (int i = 0; i < N+1; i++) {
            a[row1][i] += k * a[row2][i];
        }
    }

    // Multiplying all entries of a row by a non-zero constant
    private void multiplyRow(int row, double k) {
        for (int i = 0; i < N+1; i++) {
            a[row][i] *= k;
        }
    }

    private void gaussianElimination() {
    /* Gaussian Elimination */
        // Forward elimination
        for (int p = 0, q = 0; p < M && q < N; q++) {
            int max = p;
            for (int i = p+1; i < M; i++) {
                if (Math.abs(a[i][q]) > Math.abs(a[max][q])) {
                    max = i;
                }
            }

            if (Math.abs(a[max][q]) > EPSILON) {
                switchRow(p, max);
                for (int i = p+1; i < M; i++) {
                    addRow(i, p, -a[i][q]/a[p][q]);
                }
                p++;
            }
        }

        echelon = true;
    }

    private void gaussJordanElimination() {
    /* Gauss-Jordan Elimination */
        gaussianElimination();

        // Back subtitution
        for (int p = M-1; p > 0; p--) {
            for (int q = 0; q < N; q++) {
                if (Math.abs(a[p][q]) > EPSILON) {
                    for (int r = 0; r < p; r++) {
                        multiplyRow(p, 1/a[p][q]);
                        addRow(r, p, -a[r][q]);
                    }
                    break;
                }
            }
        }
        multiplyRow(0, 1/a[0][0]);

        reduced = true;
    }

    public void solve() {
        gaussianElimination();

        int row = 0;
        boolean isConsistent = true;

        formatZeroRow();

        show();

        loopInitialization:
        for (int p = M-1; p >= 0; p--) {
            for (int q = 0; q <= N; q++) {
                if (Math.abs(a[p][q]) > EPSILON) {
                    row = p;
                    isConsistent = q != N;
                    break loopInitialization;
                }
            }
        }

        if (isConsistent) { // Has solution
            if (row < N-1) { // Infinitely many solutions
                System.out.println("This system has infinitely many solutions");
                if (reduced) {
                    boolean[] isFreeVariable = new boolean[N];

                    for (int i = 0, j = 0; j < N; j++) {
                        if (Math.abs(a[i][j]) > EPSILON) {
                            isFreeVariable[j] = false;
                            i++;
                        } else {
                            isFreeVariable[j] = true;
                        }
                    }

                    for (int i = 0, j = 0; j < N; j++) {
                        if (isFreeVariable[j]) {
                            System.out.printf("x%d = %c", j+1, 'p'+j);
                        } else {
                            System.out.printf("x%d = %f", j+1, a[i][N]);
                            for (int k = j+1; k < N; k++) {
                                if (Math.abs(a[i][k]) > EPSILON)
                                    System.out.printf("%+f%c", -a[i][k], 'p'+k);
                            }
                            i++;
                        }
                        System.out.println();
                    }
                } else {
                    double[][] x = new double[N][N+1];
                    
                    for (int i = N-1; i >= 0; i--) {
                        for (int j = row; j <= i; j++) {
                            if (Math.abs(a[row][j]) > EPSILON) {
                                if (j < i) {
                                    x[i][i] = 1;
                                } else {
                                    x[i][N] = a[row][N];
                                    
                                    for (int k = j+1; k < N; k++) {
                                        for (int r = k; r <= N; r++) {
                                            x[i][r] -= a[row][k] * x[k][r];
                                        }
                                    }
                                    
                                    for (int k = j+1; k <= N; k++) {
                                        x[i][k] /= a[row][j];
                                    }
        
                                    row--;
                                }
                                break;
                            }
                        }
                    }
                    
                    for (int i = 0; i < N; i++) {
                        System.out.printf("x%d = % f", i+1, x[i][N]);
                        for (int j = 0; j < N; j++) {
                            if (Math.abs(x[i][j]) > EPSILON)
                                System.out.printf("%+f%c", x[i][j], 'p'+j);
                        }
                        System.out.println();
                    }
                }
            } else { // One unique solution
                System.out.println("This system has one unique solution");
                if (reduced) {
                    for (int i = 0; i < N; i++) {
                        System.out.printf("x%d = %f%n", i+1, a[i][N]);
                    }
                } else {
                    double[] x = new double[N];
                    
                    for (int i = row; i >= 0; i--) {
                        x[i] = a[i][N];
                        for (int j = i+1; j < M; j++) {
                            x[i] -= a[i][j] * x[j];
                        }
                        x[i] /= a[i][i];
                    }
        
                    for (int i = 0; i < M; i++) {
                        System.out.printf("x%d = %f%n", i+1, x[i]);
                    }
                }
            }
        } else { // No solution
            System.out.println("This system has no solution");
        }
    }

    public double interpolate(double x) {
        if (!reduced) {
            gaussJordanElimination();
        }

        double f = 0;
        for (int i = 0; i < M; i++) {
            f += Math.pow(x, i) * a[i][N];
            System.out.printf("a%d = %f%n", i, a[i][N]);
        }
        
        return f;
    }

    public void show() {
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {
                System.out.printf("%8.3f ", a[i][j]);
            }
            System.out.printf("| %8.3f%n", a[i][N]);
        }
        System.out.println();
    }

    private void formatZeroRow() {
        loopFormatZeroRow:
        for (int p = M-1; p >= 0; p--) {
            boolean zeroRow = true;
            for (int q = N-1; q >= p; q--) {
                if (Math.abs(a[p][q]) > EPSILON) {
                    zeroRow = false;
                    break;
                }
            }
            if (zeroRow) {
                if (Math.abs(a[p][N]) > EPSILON) {
                    for (int i = p-1; i >= 0; i--) {
                        for (int j = N-1; j >= 0; j--) {
                            if (Math.abs(a[i][j]) > EPSILON) {
                                a[i+1][N] = 1;
                                for (int k = i+2; k < M; k++) {
                                    a[k][N] = 0;
                                }
                                break loopFormatZeroRow;
                            }
                        }
                    }
                }
            } else {
                break;
            }
        }
    }

    public static void test(double[][] A, double[] b) {
        SystemOfLinearEquation gaussian = new SystemOfLinearEquation(A, b);
        gaussian.solve();
    }

    // 3-by-3 nonsingular system
    public static void test1() {
        double[][] A = {
            { 0, 1,  1 },
            { 2, 4, -2 },
            { 0, 3, 15 }
        };
        double[] b = { 4, 2, 36 };
        test(A, b);
    }

    // 3-by-3 nonsingular system
    public static void test2() {
        double[][] A = {
            {  1, -3,   1 },
            {  2, -8,   8 },
            { -6,  3, -15 }
        };
        double[] b = { 4, -2, 9 };
        test(A, b);
    }

    // 5-by-5 singular: no solutions
    // y = [ -1, 0, 1, 1, 0 ]
    public static void test3() {
        double[][] A = {
            {  2, -3, -1,  2,  3 },
            {  4, -4, -1,  4, 11 },
            {  2, -5, -2,  2, -1 },
            {  0,  2,  1,  0,  4 },
            { -4,  6,  0,  0,  7 },
        };
        double[] b = { 4, 4, 9, -6, 5 };
        test(A, b);
    }

    // 5-by-5 singular: infinitely many solutions
    public static void test4() {
        double[][] A = {
            {  2, -3, -1,  2,  3 },
            {  4, -4, -1,  4, 11 },
            {  2, -5, -2,  2, -1 },
            {  0,  2,  1,  0,  4 },
            { -4,  6,  0,  0,  7 },
        };
        double[] b = { 4, 4, 9, -5, 5 };
        test(A, b);
    }

    // 3-by-3 singular: no solutions
    // y = [ 1, 0, 1/3 ]
    public static void test5() {
        double[][] A = {
            {  2, -1,  1 },
            {  3,  2, -4 },
            { -6,  3, -3 },
        };
        double[] b = { 1, 4, 2 };
        test(A, b);
    }

    // 3-by-3 singular: infinitely many solutions
    public static void test6() {
        double[][] A = {
            {  1, -1,  2 },
            {  4,  4, -2 },
            { -2,  2, -4 },
        };
        double[] b = { -3, 1, 6 };
        test(A, b);
    }

    public static void test7() {
        double[][] A = {
            {  1,  2,  3,  4,  5 },
            {  0,  0,  5,  4,  3 },
            {  0,  0,  0,  0,  0 },
            {  0,  0,  5,  4,  3 },
        };
        double[] b = { -3, 1, 0, 3 };
        test(A, b);
    }

    public static void example1() {
        double[][] A = {
            { 0.31, 0.14, 0.30, 0.27 }, 
            { 0.26, 0.32, 0.18, 0.24 }, 
            { 0.61, 0.22, 0.20, 0.31 }, 
            { 0.40, 0.34, 0.36, 0.17 },
        };
        double[] b = {1.02, 1.00, 1.34, 1.27};
        SystemOfLinearEquation gaussian = new SystemOfLinearEquation(A, b);
        gaussian.solve();
    }

    public static void example2() {
        double[][] A = {
            {   1,   7,  -2,   0,   8 }, 
            {   1,   7,  -1,   4,   0 }, 
            {   2,  14,  -4,   1, -13 }, 
            {   2,  14,  -4,   0,  16 }
        };
        double[] b = {-3, 2, 3, -6};
        SystemOfLinearEquation gaussian = new SystemOfLinearEquation(A, b);
        gaussian.solve();
    }

    public static void example3(int n) {
        double[][] H = hilbert(n);
        double[] b = ones(n);
        SystemOfLinearEquation gaussian = new SystemOfLinearEquation(H, b);
        gaussian.solve();
    }

    public static void example4() {

    }

    public static void example5() {
        
    }

    public static void example6(double a, double b, int n) {
        double[] xs = new double[n+1];
        double[] fs = new double[n+1];
        double h = (b-a)/n;
        for (int i = 0; i <= n; i++) {
            xs[i] = a+i*h;
            fs[i] = Math.exp(-xs[i])/(1+Math.pow(xs[i], 1/2)+Math.pow(xs[i], 2));
        }

        SystemOfLinearEquation gaussian = new SystemOfLinearEquation(xs, fs);

        gaussian.gaussJordanElimination();
        
        System.out.printf("p(x) = %f", gaussian.a[0][gaussian.N]);
        for (int i = 1; i <= n; i++) {
            System.out.printf("%+fx^%d", gaussian.a[i][gaussian.N], i);
        }
        System.out.println();
    }

    public static void example7(double x) {
        double[] xs = {0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3};
        double[] fs = {0.003, 0.067, 0.148, 0.248, 0.370, 0.518, 0.697};

        SystemOfLinearEquation gaussian = new SystemOfLinearEquation(xs, fs);

        System.out.printf("f(%f) = %f%n", x, gaussian.interpolate(x));
    }

    public static void example8(int tahun) {
        double[] xs = {1950, 1955, 1960, 1965, 1966, 1967, 1968, 1969};
        double[] fs = {33525, 46519, 53941, 72319, 75160, 76160, 84690, 90866};

        SystemOfLinearEquation gaussian = new SystemOfLinearEquation(xs, fs);
        
        System.out.printf("Prediksi harga rumah baru pada tahun %d $%f juta%n", tahun, gaussian.interpolate(tahun));
    }

    public static void example9(double t) {
        double[] xs = {40, 50, 60, 70, 80, 90};
        double[] fs = {1.66, 1.41, 1.22, 1.06, 0.93, 0.84};

        SystemOfLinearEquation gaussian = new SystemOfLinearEquation(xs, fs);
        
        System.out.printf("Taksiran viskositas pada suhu %fderajat Fahrenheit %fft^2/detik%n", t, gaussian.interpolate(t));
    }

    public static double[][] hilbert(int n) {
        double[][] H = new double[n][n];
        for (int i = 0; i < n; i++) 
            for (int j = 0; j < n; j++) 
                H[i][j] = 1/(i+j-1);
        return H;
    }

    public static double[] ones(int n) {
        double[] one = new double[n];
        for (int i = 0; i < n; i++) 
            one[i] = 1;
        return one;
    }

    // sample client
    public static void main(String[] args) {

        try                { test1();             }
        catch(Exception e) { e.printStackTrace(); }
        System.out.println("--------------------------------");

        try                { test2();             }
        catch(Exception e) { e.printStackTrace(); }
        System.out.println("--------------------------------");

        try                { test3();             }
        catch(Exception e) { e.printStackTrace(); }
        System.out.println("--------------------------------");

        try                { test4();             }
        catch(Exception e) { e.printStackTrace(); }
        System.out.println("--------------------------------");

        try                { test5();             }
        catch(Exception e) { e.printStackTrace(); }
        System.out.println("--------------------------------");

        try                { test6();             }
        catch(Exception e) { e.printStackTrace(); }
        System.out.println("--------------------------------");

        try                { test7();             }
        catch(Exception e) { e.printStackTrace(); }
        System.out.println("--------------------------------");
    }

}