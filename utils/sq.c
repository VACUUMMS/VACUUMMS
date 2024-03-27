/* sq.c */

#include <vacuumms/std.h>
#include <vacuumms/rng2.h>
#include <vacuumms/param.h>
#include <errno.h>

// In:  x, out: x^2 

int main(int argc, char* argv[])
{
  double x, y, z, d;
  char line[80];
  char *xs;

  while (TRUE)
  {
    fgets(line, 80, stdin);
    if (feof(stdin)) break;

    xs = strtok(line, "\n");

    x = strtod(xs, NULL);

    printf("%lf\n", x*x);
  }
}

