# VACUUMMS variational module

The variational module is an extension to VACUUMMS that allows for finding the paths of least energy along the insertion energy landscape. 

<image of general problem>

Taking the set of pairs of Voronoi vertices (in lieu of CESA cavity centers, which are nearly identical but which lack the implicit pairing of Voronoi edges) as input. The following protocol is applied for each edge:

- End points and an energy function (or molecular configuration) are taken as input. A set of evenly spaced points (the variational points) is created between these endpoints, comprising the variational path. 

- An iterative gradient descent is performed for each of the variational points, nudging each point to a lower energy. This nudging is perpendicular to the variational path, and is determined by evaluating directional derivatives of the energy in the plane (or line in 2D) perpendicular to the path:

  <image of perpendicular>

  The direction of these derivatives is derived by rotating the (x, y, z) system such that the z-axis points in the direction of the variational path at the variational point being evaluated. The parameters of the rotation are determined from mapping z, such that the same rotation maps unit vectors pointing along the x and y axes to mutually orthogonal vectors, orthogonal to the variational path at that point. The **gradient** of the energy function in the plane spanned by the resulting vectors is calculated and used to generate the perturbation of the variational point. (see section on quaternion info)

- Calling the `iterate()` function causes this process is repeated for each of the variational points along the curve, each being evaluated individually before updating the entire curve to the perturbed points. 

- Since these perturbations invariably cause an uneven stretching of the curve, a `rebalancePoints()` method is provided that performs an adaptive re-spacing of points evenly along the variational curve:

  <image of respacing>

- Convergence criteria have not yet been implemented. A likely candidate is when mean-square-displacement of variational points falls below a threshold value with each successive iteration, or alternatively, a fixed number of iterations is performed. 

## Quaternion form for gradient descent in perpendicular plane

Each variational point *p* is neighbored by a point before and a point after (one of being the start or end point of the variational curve itself for the first and last variational points, respectively). These two neighboring points are used to define an axis perpendicular (within numerical approximation) to the variational curve at the variational point of interest. Let these points be represented as **p_fore** and **p_aft** and the unit vector tangent to the variational curve as **n**:

​		**n = (p_aft - p_fore)** /  **| p_aft - p_fore |**  

​		<image of this>

as the axis about which the base coordinate system is rotated through an angle *theta* resulting in the unit vector **k** mapping to **n**. In this way, the **i** and **j** unit vectors map to unit vectors orthogonal to the variational curve, as well as orthogonal to one another, thereby providing a basis set that spans the plane perpendicular to the variational curve. In this plane, a step of gradient descent is performed for each variational point. Only after all of the perturbations have been calculated for each individual point is the entire curve updated. 

The cross product **k** X **n** describes the rotation from **k** to **n** through an angle *theta* around the axis **u**:

​		**k** X **n** = |**k**||**n**|sin(*theta*)

										 |	 i	 j	 k	 |
			k X n = |k||n|sin(theta) =	 |	 0	 0	 1	 |
										 |	n_x	n_y	n_z	 | 

**u** is the unit vector pointing in the direction of the cross product:

​		**u** = **k** X **n** / |**k** X **n**|

and *theta* is extracted *via* the arcsine function, with the dot-product **k** * **n** preserving the sense (+/-) of rotation:

​		*theta* = *arcsin*[ ( |**k** X **n**| /|**k**||**n**| ) ( **k** * **n** /|**k** * **n**| ) ] 

Thus the quaternion-based mapping of the rotation is:

​		***q** = cos(theta/2) + **u** sin(theta/2)*

​		**v'** = **qvq*** 

where **v** is the original vector and **v'** is the vector in the rotated coordinate system, and **q*** is the quaternion conjugate of **q**. Thus, **k** rotates to **n** and by definition, **i** and **j** rotate to **directional_x** and **directional_y** respectively. These directionals are used to evaluate the directional gradient in the potential field E and thus generate perturbations:

​		**p_n+1** = **p_n** + *alpha* * ***grad*** E(**p_n**)

### Example 1

Consider a variational point where the variational curve points in the direction (1,1,1). Let *a* = *sqrt(*3*)*/3 and b = sqrt(2). The unit vector **n** is thus (*a, a, a*). Mapping **v** = **k** (unit vector in z-direction) to **v' = n** means evaluating the cross product determinant:

					  |	 i	 j	 k	 |
			k X n =   |	 0	 0	 1	 | = (-a, a, 0)
					  |	 a	 a	 a	 | 

The magnitude of which is:

​		|**k** X **n**| = *sqrt( -a* * -*a* + *a* * *a* ) = *a sqrt*(2) = *ab*			

Since the dot product **k** * **n** is *a*, which is positive, the value for *theta* can be taken directly from the arcsine (dagger: there is an additional step if the dot product is negative):

​		*theta* = *arcsin ( ab )* = 0.955.

and the rotation vector **u** is:

​		**u** = k X n / |k X n| = (*-a, a, 0*) */ab* = (-1, 1, 0) / *b*.

Thus giving the rotation quaternion conjugate pair:

​		**q** = *cos*(*theta/2*) + **u** sin(*theta/2*) = 0.888 - 0.325**i** + 0.325**j**

​		**q*** = *cos*(*theta/2*) - **u** sin(*theta/2*) = 0.888 + 0.325**i** - 0.325**j**

Rotating the **i** and **j** vectors gives directionals:

​		**directional_x** =  **qiq*** = 0.788**i** + 0.212**j** + 0.577**k**

​		**directional_y** = **qjq*** = -0.212**i** + 0.788**j** - 0.577**k**

And as a sanity check, rotating the unit vector **k** should return the original vector **n**, within the limits of numerical precision. Preserving three significant figures here returns:

​		**qkq*** = 0.577**i** + 0.577**j** + 0.576**k** ~ (*a, a, a*)

### Example 2

Consider a variational point where the variational curve points in the direction **n** = (0,1,0) (parallel to y-axis). Mapping **v** = **k** (unit vector in z-direction) to **v'** = **n** means evaluating the cross product determinant:

					|	i	j	k	|
		k X n = 	|	0	0	1	|	=	(-1, 0, 0)
					|	0	1	0	| 

The magnitude of which is:

		|k X n| = 1				

thus:

		theta = arcsin(k X n) = arcsin(-1) = - pi/2

and:

		u = k X n / |k X n| = (-1, 0, 0) = -i

Thus giving the rotation quaternion conjugate pair:

		q = cos(pi/4) - i sin(pi/4) = (sqrt(2)/2) (1 - i)
		q* = cos(pi/4) + i sin(pi/4) = (sqrt(2)/2) (1 + i)

Applying this to the vector **v** = (1, 0, 0):

```	**v'** = **qvq*** = (sqrt(2)/2) * (sqrt(2)/2) (1 - i) (i) (1 + i)
	v' = qvq* = (sqrt(2)/2) * (sqrt(2)/2) (1 - i) (i) (1 + i)
	          = (1/2) (i - (-1))(1 + i)
              = (1/2) (i + (-1) + 1 + i)
              = i
```



# Prerequisites

The variational module requires the voronoi module (which in turn requires voro++) to be installed. This integration is handled by CMake and/or spack, but does require that voro++ be available. 

# Developer overview

The variational module is written in C++, and makes use of existing VACUUMMS code via *'extern "C" '* declarations. Relevant classes are Variational2D and Variational3D which take a set of endpoints, number of variational points to be used, and either a molecular configuration or alternatively, a pointer to a function of x, y, (and z in 3D case) and returns a scalar energy value.


```
class Variational2D
{
    private:

        float start_x;
        float start_y;
        float end_x;
        float end_y;
        int n_var_points;
        Configuration *configuration;
        float(*energy_function)(float x, float y);

        bool use_configuration_energy;

        float* var_x;
        float* var_y;

    public:

        void init(float start_x,
             float start_y,
             float end_x,
             float end_y,
             int n_var_points);

        Variational2D(float start_x,
                      float start_y,
                      float end_x,
                      float end_y,
                      int n_var_points,
                      float(*energy_function)(float x, float y));

        Variational2D(float start_x,
                      float start_y,
                      float end_x,
                      float end_y,
                      int n_var_points,
                      Configuration *c);

        void printValues();
        void iterate();
        void iterateWork();

        ~Variational2D();

        float rebalancePoints2D();

};

class Variational3D
{
    private:

        float start_x;
        float start_y;
        float start_z;
        float end_x;
        float end_y;
        float end_z;
        int n_var_points;
        Configuration *configuration;
        float(*energy_function)(float x, float y, float z);

        bool use_configuration_energy;

        float* var_x;
        float* var_y;
        float* var_z;

    public:

        void init(float start_x,
             float start_y,
             float start_z,
             float end_x,
             float end_y,
             float end_z,
             int n_var_points);

        Variational3D(float start_x,
                      float start_y,
                      float start_z,
                      float end_x,
                      float end_y,
                      float end_z,
                      int n_var_points,
                      float(*energy_function)(float x, float y, float z));

        Variational3D(float start_x,
                      float start_y,
                      float start_z,
                      float end_x,
                      float end_y,
                      float end_z,
                      int n_var_points,
                      Configuration *c);

        void printValues();
        void iterate();
        void iterateWork();

        ~Variational3D();

        float rebalancePoints3D();

};
```

