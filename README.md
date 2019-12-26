## 2018 June

1. Let A = [1:5; 6:10; 16:20] be a matrix in **MATLAB**. Write a <u>single</u> **MATLAB** command for each of the following:

   (a)	Reverse the order of even columns of A.<span style="float:right">1 mark</span>

   ```matlab
   A(:, 2:2:end) = A(:, 4:-2:1)
   ```

   (b)	Find the sum of all entries of A in even columns.<span style="float:right">2 marks</span>

   ```Matlab
   sum(sum(A(:, 2:2:end)))
   ```

   (c)	Insert 11:15 as a new third row of A.<span style="float:right">2 marks</span>

   ```matlab
   A = [A(1:2, :); 11:15; A(3:end, :)]
   ```

   (d)	Remove the second and the second last columns of A.<span style="float:right">2 marks</span>

   ```matlab
   A(:, [2 end-1]) = []
   ```

   (e)	Increase each entry of A which is less than 20 and greater than 3 by 3.<span style="float:right">2 marks</span>
   
   ```matlab
   A(A < 20 & A > 3) = A(A < 20 & A > 3) + 3
   ```
   
2. (a)	Consider the sequence of integers $\{F_n\}_{n=1}^{∞}$ where F~1~ = 1, F~2~ = 3sin(45°), F~3~=7cos(120°) and 
   $$
   F_n = F_{n-1} - 2F_{n-2} + F_{n-3}
   $$
   for n≥4.

   ​	(i)	Write a function m-file called F.m such that given any positive integer k, F(k) returns the row vector [F~1~ F~2~ ... F~k~].<span style="float:right">3 marks</span>

   ```matlab
function out=F(k)
     out = [1 3*sin(pi*45/180) 7*cos(pi*120/180)];
    for i=4:k
       out(end+1) = out(end) - 2*out(end-1) + out(end-2);
     end
     out = out(1:k);
   end
   ```
   
   ​	(ii)	Write a script to determine the largest positive integer n such that $F_1+F_2+...+F_n<10000$<span style="float:right">3 marks</span>
   
   ```matlab
   n = 1;
   while (sum(F(n))) < 10000
     n = n+1;
   end
   n = n-1;
   ```
   
   (b)	Calculate the value of $f_{120}$ using the <u>recurrence relation</u>(递推关系式), $f_1=1, f_2=\frac{1}{4}, f_n=f_{n-1}-\frac{3}{4}, n=3,4,...$<span style="float:right">3 marks</span>
   
   ```matlab
   f = [1];
   for i=2:120
     f(end+1) = f(end) - 3/4;
   end
   f120 = f(120);
   ```
   
3. (a)	Given the functions $u(x) = \frac{sin(x)}{x}, v(x)=\frac{sin(x)}{x^2} - \frac{cos(x)}{x}$ for $x=0.001, 0.01, 0.02, ..., 5.0$, and the function $w(x) = 3\frac{v(x)}{x} - u(x)$.<span style="float:right">10 marks</span>

   ​	(i)	Write a script to generate the graphs of u(x), v(x) and w(x) on the same figure with a stepsize of 0.01, Labels the figures completely.

   ​	(ii)	Use the function <u>fzero</u> to solve the equation $v(x) = u(x)$ and $w(x)=u(x)$. Compare the values of these roots with the graph in (i).

   ```matlab
   % question 01
   x = 0:0.01:5;
   u = @(x) sin(x) ./ x;
   v = @(x) sin(x) ./ (x.^2) - cos(x) ./ x;
   w = @(x) 3*v(x)./x - u(x);
   plot(x, [u(x)' v(x)' w(x)']);
   xlabel('x'); ylabel('y');
   legend('u(x)', 'v(x)', 'w(x)');
   % question 02
   fun1 = @(x) u(x)-v(x);
   fun2 = @(x) w(x)-u(x);
   root1 = fzero(fun1, 2)
   root2 = fzero(fun2, 2.5)
   ```

   (b)	Write a script using the subplot command to plot
   $$
   y_1(x) = e^{-2x}+3sin(0.5x) \\
   y_2(x) = \frac{d}{dx}(xe^{-0.25x})
   $$
   ​		for the interval [0, 2π]. You may use a stepsize of 0.1, and label the figures completely.<span style="float: right">8 marks</span>

   ```matlab
   x = 0:0.1:2*pi;
   y1 = @(x) exp(-2*x) + 3*sin(0.5*x);
   y2 = @(x) exp(-0.25*x) - 0.25*x.*exp(-0.25*x);
   
   ```

   (c)	Consider the function 
   $$
   r(t) = e^{cos(t)} - 2cos(4t) + sin(\frac{t^5}{12})
   $$
   ​		over the domain 0 ≤ t ≤ 2π.

   ​		(i)	Plot r(t) in the polar coordinate system.

   ​		(ii)	By using a suitable MATLAB command, obtain the cartesian 				coordinates (x, y) and plot y against x.

4. (a)	Define the polynomial 
   $$
   P(x) = 6x^5 + 3x^3 - 15x^2 - 2x + 9
   $$
   ​		and perfome the following:

   ​	(i)	Evaluate $\frac{d^3P(x)}{dx^3}$

   ​	(ii)	Solve $P(x)=0, \frac{d^3P(x)}{dx^3}=0$

   ​	(iii)	Evaluate $(P(x))^3$

   (b)	Write **MATLAB** statements to calculate the doubles sum
   $$
   S = \sum^{20}_{n=1}\sum^{30}_{m=1}(n^2+(-2)^nm^2)
   $$


## 2019 January

1. (a)	if u = [2, 4, 6, ..., 120] and v = [1, 3, 5, ..., 120] are two vectors, then write the **MATLAB** commands to find the following

   ​	(i)	Represent the vectors u and v using colon notation.

```matlab
u = [2:2:120];
v = [1:2:120];
```

​			(ii)	Length of the vectors u and v.

```matlab
fprintf("the length of u is %d", length(u));
fprintf("the length of u is %d", length(v));
```

​			(iii)	Find the last element of the two vectors.

```matlab
fprintf("the last element of u is %d", u(end));
fprintf("the last element of v is %d", v(end));
```

​			(iv)	Find the inner product of the vectors u and v.

```matlab
fprintf("the inner product of u and v is %d", sum(u .* v));
```

​			(v)	Find the square of the vectors u and v using Hadamard Product.

```matlab
u .* v;
```

​			(vi)	Find the division of the vector u by u and also division of the vector v by v.

```matlab
u ./ u;
v ./ v;
```

​	(b)	Consider the following functions:
$$
y_1 = 2cosx \\
y_2 = 2sinx \\
y_3 = 1 +cosx \\
y_4 = 1 + sinx
$$

​			 (i)	Write **MATLAB** script file to plot all the functions in the same axes in the interval 0≤x≤2π.

```matlab
x = 0:0.01:2*pi;
y1 = 2*cos(x);
y2 = 2*sin(x);
y3 = 1+cos(x);
y4 = 1+sin(x);
plot(x, [y1' y2' y3' y4']);
```

​			 (ii)	Write **MATLAB** statements to get subplot and label completely.

```Matlab
x = 0:0.01:2*pi;
y1 = 2*cos(x);
y2 = 2*sin(x);
y3 = 1+cos(x);
y4 = 1+sin(x);
y = [y1' y2' y3' y4'];
for i=1:4
    subplot(4, 1, i);
    plot(x, y(:, i));
    xlabel('x');
    ylabel('y');
end
```

2. (a)	Given $A =\begin{pmatrix} 1&2&3 \\ 1&5&6 \\ 7&8&9 \end{pmatrix}$ and $B =\begin{pmatrix} 2&2&3&1 \\ 2&5&6&2 \\ 6&3&9&3 \end{pmatrix}$ </br>Write down the **MATLAB** statements to find the following:

   ​		(i)	Size of the mitrices A and B.

```matlab
size_a = size(A);
size_b = size(B);
```

​			(ii)	Convert the matrices A to identity matrix.

```matlab
A = eye(size_a);
```

​			(iii)	Convert the matrices A to diagonal matrix(对角矩阵).

```matlab
A = diag(diag(A));
```

​			(iv)	Product of two matrices A and B.

```matlab
A * B
```

​			(v)	Determinant of the matrices(矩阵的行列式) A.

```matlab
det(A)
```

​			(vi)	Inverse of the matrices(逆矩阵) A.

```matlab
inv(A)
```

​			(vii)	Extract a matrix containing second column of the matrix A and second column of matrix B.

```matlab
[A(:, 2) B(:, 2)]
```

​	(b)	Write down the **MATLAB** statement to solve the following system of linear equations using Cramer's rule, 
$$
x + 2y + 3z = 1\\
3x + 3y + 4z = 1 \\
2x + 3y + 3z = 2
$$
where the Cramer's rule is given by: 
$$
\text{x} = \frac{|A_x|}{|A|}, \text{y} = \frac{|A_y|}{|A|}, \text{z} = \frac{|A_z|}{|A|},
$$

```matlab
k = [1 2 3; 3 3 4; 2 3 3];
y = [1; 1; 2];
x = det([y k(:, 2:3)]) / det(k);
y = det([y k(:, 1) k(:, 3)]) / det(k);
z = det([y k(:, 1:2)]) / det(k);
```

​	(c)	Write the **MATLAB** function program to find the area of triangle whose sildes are 10, 15 and 20 unit. The area of triangle is given by: 
$$
Area = \sqrt{s(s-a)(s-b)(s-c)}, \text{and } s = \frac{a+b+c}{2}
$$

```matlab
function out=area(a, b, c)
  s = (a+b+c) / 2;
  out = sqrt(s*(s-a)*(s-b)*(s-c));
end
```

3. (a)	==TODO:== Consider the following equations,
   $$
   y^2 = 10 -x \text{ and } x = (y-2)^2
   $$
   ​	(i)	Plot the equations using **MATLAB** anoymous function for -10≤x≤20 and -2≤y≤4.

```matlab
y^2 = @(x) (10-x);
x = -10:0.1:20;
plot(x, y^2);
```

   ​	(ii)	Hence, determine the intersection point(s) between the equations.

   (b)	==TODO:== Determine

   ​	(i)	the area of the region bounded by $y^2 =10-x$ and $x=(y-2)^2$

   ​	(ii)	the volume of the region bounded by $y^2 =10-x$ and $x=(y-2)^2$ when they are revolved 2π about y-axis.

   Give your answer in rational numbers.

4. (a)	Find a fourth degree polynomial that passes through (-1, 10), (1, 4), (2, 10), (3, 22), and (4, 10). Hence, plot the polynomial.

```matlab 
% solution1 - use solver, cheat!
syms a b c d e
f = @(x) a*x^4 + b*x^3 + c*x^2 + d*x + e
eqns = [f(-1) == 10, f(1) == 4, f(2) == 10, f(3) == 22, f(4) == 10];
vars = [a b c d e];
[sol_a, sol_b, sol_c, sol_d, sol_e] = solve(eqns,vars)
% a=-1, b=5, c=-2, d=-8, e=10
fun = @(x) (-1)*x.^4 + 5*x.^3 + (-2)*x.^2 + (-8)*x + 10;
fplot(fun, [-1, 4]);

% solution2 - use polyfit
x = [-1 1 2 3 4];
y = [10 4 10 22 10];
p = polyfit(x, y, 4);
x1 = linspace(-1, 4, 5000);
y1 = polyval(p, x1);
plot(x, y, 'o');
hold on;
plot(x1, y1);
hold off;
```

   (b)	Find the local maximum and local minimum points for the polynomial in part (a) for the range of -1≤x≤4 accurate to four decimal places.

```matlab
fun = @(x) (-1)*x.^4 + 5*x.^3 + (-2)*x.^2 + (-8)*x + 10;
fplot(f, [-1, 4]);
x = -1:0.001:4;
xmin = fminbnd(fun, -1, 4)
xmax = x(max(f(x)) == f(x))
```

---

## 2019 June

1. (a)	Write a **MATLAB** script M-file to find the roots of the quadratic equation $ax^2+bx+c=0$. Then, find the roots of 

$$
x^2 + 4x+4=0 \tag{i}
$$
$$
x^2+5x+6=0 \tag{ii}
$$
$$
x^2+x+1=0 \tag{iii}
$$

```matlab
out1 = quad_roots(1, 4, 4);
out2 = quad_roots(1, 5, 6);
out3 = quad_roots(1, 1, 1);

function out=quad_roots(a, b, c)
  out = roots([a b c]);
end
```

​		(b)	Write **MATLAB** statements to calculate the sum(S~1~)+sum(S~2~) where 
$$
S_1=1:2:100 \quad and \quad S_2=-2:-2:-100
$$

```matlab
s1 = 1:2:100;
s2 = -2:-2:-100;
sum(sum(s1) + sum(s2))
```


2. (a)	Let $y_1 = e^{-1.5x}sin(10x)$ and $y_2 = e^{1.5x}cos(10x); 0\leq x \leq 2\pi$

   ​		(i)	Plot y~1~ and y~2~ in the same figure.

   ​		(ii)	plot y~1~ and y~2~ in the separate figure.

   ​		Give the labels and legends for the two figures.

```matlab 
y1 = @(x) exp(-1.5*x).*sin(10*x);
y2 = @(x) exp(1.5*x).*cos(10*x);
fplot(y1, [0, 2*pi]);
hold on;
fplot(y2, [0, 2*pi]);
hold off;
xlabel('x');
ylabel('y');
legend('y_1=e^{-1.5*x}sin(10x)', 'y_2=e^{1.5*x}cos(10x)');
```

   (b)	Write down **MATLAB** statement to draw contour plot, meshplot, surfaceplot $sin(xy^2)+cos(yx^2)$ with $x$ in the range $[-2\pi, 2\pi]$ and $y$ in the range $[0, 4\pi]$. Give the title and labels.

```matlab
[x, y] = meshgrid([-2*pi, 2*pi], [0, 4*pi]);
z = sin(x.*y.^2) + cos(y.*x.^2);
contour(x, y, z);
title('contour plot of sin(xy^2)+cos(yx^2)');
xlabel('x'); ylabel('y'); zlabel('z');
mesh(x, y, z);
surf(x, y, z);
```

3. (a) 	Write a **MATLAB** program that reads the user's age and then print "You are a child" if the age < 18, "You are an adult." if 18 ≤ age ≤ 65, and "You are a senior citizen" if age ≥ 65.

```matlab
function out=read_age(age)
  if age<18
    printf('You are a child.');
  elseif age>=18 && age<=65
    printf('You are an adult.');
  else
    printf('You are a senior citizen.');
  end
end
```

   (b)	Write a **MATLAB** function to find Fibonacci sequence F(1) = 0, F(2) = 1, F(n) = F(n-1) + F(n-2) using loop structure. 

```matlab
function out=fib(n)
  out=[0 1];
  for i=3:n
    out(end+1)=out(end)+out(end-1);
  end
  out = out(1:n);
end
```

4. (a)	Write a **MATLAB** function to find the maximum and minimum of given three numbers a, b ,c.

```matlab
function [maximum, minimum] = max_and_min(a, b, c)
  maximum = max([a b c]);
  minimum = min([a b c]);
end
```

   (b)	==TODO:== Given $f(x) = x^5 - 3x^4 + 3x^3 - 2x^2 - 5.$ Write the **MATLAB** commands to find 

   ​		(i)	roots of $f(x)$

   ​		(ii)	derivative of $f(x)$ 

   ​		(iii)	min and max of $f(x)$

   ​		(iv)	critival points of $f(x)$

```matlab
% question 01
r = roots([1 -3 3 -2 -5]);
% question 02
syms x
f = x^5 - 3*x^4 + 3*x^3 - 2*x^2 - 5;
diff(f)
% question 03
```

5. (a)	Determine a polynomial that passes through the points (-1, 1), (0, -2), (1, 0) and (2, -1) by using **MATLAB**.

```matlab
x = [-1 0 1 2];
y = [1 -2 0 -1];
p = polyfit(x, y, 3);
x1 = linspace(-1, 2, 1000);
y1 = polyval(p, x1);
plot(x, y, 'o');
hold on;
plot(x1, y1);
```

   (b)	Find the area between the curves $f_1(y) = 3+y^2$, $f_2(y)=2-y^2$, $y=1$, $y=-2$ by using **MATLAB**.

```matlab
y = linspace(-4, 3, 1000);
f1 = 3 + y.^2;
f2 = 2 - y.^2;
plot(y, [f1', f2']);
hold on;
xline(-2);
xline(1);
```