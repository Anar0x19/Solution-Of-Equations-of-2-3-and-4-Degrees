// This program can solve quadratic equations, cubic equations, and 4th degree equations

#include <stdio.h>
#include <math.h>
#include <complex.h>

void solve_quadratic(double a, double b, double c, double complex roots[2]) {
    double discriminant = b * b - 4 * a * c;
    
    if (discriminant >= 0) {
        roots[0] = (-b + sqrt(discriminant)) / (2 * a);
        roots[1] = (-b - sqrt(discriminant)) / (2 * a);
    } else {
        roots[0] = (-b + csqrt(discriminant)) / (2 * a);
        roots[1] = (-b - csqrt(discriminant)) / (2 * a);
    }
}

void solve_cubic(double a, double b, double c, double d, double complex roots[3]) {
    b /= a;
    c /= a;
    d /= a;

    double delta_0 = b * b - 3 * c;
    double delta_1 = 2 * b * b * b - 9 * b * c + 27 * d;

    double complex C = cpow((delta_1 + csqrt(delta_1 * delta_1 - 4 * delta_0 * delta_0 * delta_0)) / 2.0, 1.0 / 3.0);

    if (cabs(C) == 0) {
        C = cpow((delta_1 - csqrt(delta_1 * delta_1 - 4 * delta_0 * delta_0 * delta_0)) / 2.0, 1.0 / 3.0);
    }

    double complex u[3];
    u[0] = 1;
    u[1] = (-1 + csqrt(-3)) / 2.0;
    u[2] = (-1 - csqrt(-3)) / 2.0;

    for (int i = 0; i < 3; ++i) {
        roots[i] = (-1.0 / 3.0) * (b + u[i] * C + delta_0 / (u[i] * C));
    }
}

void solve_quartic(double a, double b, double c, double d, double e, double complex roots[4]) {
    b /= a;
    c /= a;
    d /= a;
    e /= a;

    double p = c - 3 * b * b / 8;
    double q = b * b * b / 8 - b * c / 2 + d;
    double r = -3 * b * b * b * b / 256 + b * b * c / 16 - b * d / 4 + e;

    double complex y_roots[3];
    solve_cubic(1, -p / 2, -r, r * p / 2 - q * q / 8, y_roots);

    double complex y = y_roots[0];
    double complex R = csqrt(b * b / 4 - y + p);
    double complex D = csqrt(3 * b * b / 4 - R * R - 2 * y + (4 * b * y - 8 * r) / (4 * R));
    double complex E = csqrt(3 * b * b / 4 - R * R - 2 * y - (4 * b * y - 8 * r) / (4 * R));

    roots[0] = (-b / 4) + (R + D) / 2;
    roots[1] = (-b / 4) + (R - D) / 2;
    roots[2] = (-b / 4) - (R + E) / 2;
    roots[3] = (-b / 4) - (R - E) / 2;
}

int main() {
    int degree;
    printf("Enter the degree of the equation (2 for quadratic, 3 for cubic, 4 for quartic): ");
    scanf("%d", &degree);

    if (degree == 2) {
        double a, b, c;
        printf("Enter coefficients a, b, c of the quadratic equation (ax^2 + bx + c = 0):\n");
        scanf("%lf %lf %lf", &a, &b, &c);

        double complex roots[2];
        solve_quadratic(a, b, c, roots);

        printf("The roots are:\n");
        for (int i = 0; i < 2; ++i) {
            printf("Root %d: %.2f + %.2fi\n", i + 1, creal(roots[i]), cimag(roots[i]));
        }

    } else if (degree == 3) {
        double a, b, c, d;
        printf("Enter coefficients a, b, c, d of the cubic equation (ax^3 + bx^2 + cx + d = 0):\n");
        scanf("%lf %lf %lf %lf", &a, &b, &c, &d);

        double complex roots[3];
        solve_cubic(a, b, c, d, roots);

        printf("The roots are:\n");
        for (int i = 0; i < 3; ++i) {
            printf("Root %d: %.2f + %.2fi\n", i + 1, creal(roots[i]), cimag(roots[i]));
        }

    } else if (degree == 4) {
        double a, b, c, d, e;
        printf("Enter coefficients a, b, c, d, e of the quartic equation (ax^4 + bx^3 + cx^2 + dx + e = 0):\n");
        scanf("%lf %lf %lf %lf %lf", &a, &b, &c, &d, &e);

        double complex roots[4];
        solve_quartic(a, b, c, d, e, roots);

        printf("The roots are:\n");
        for (int i = 0; i < 4; ++i) {
            printf("Root %d: %.2f + %.2fi\n", i + 1, creal(roots[i]), cimag(roots[i]));
        }
    } else {
        printf("Invalid degree entered. Please enter 2 for quadratic, 3 for cubic, or 4 for quartic equations.\n");
    }

    return 0;
}

