# LinearAlgebra
Theorems and Definitions from 1.1, 1.2, 1.3, ~~1.4, 1.5, 1.7, 1.8, 1.9, 2.1, 2.2,~~ 2.3, 2.8, 2.9, 3.1, 3.2, 5.1, 5.2, 5.3, 4.1, 4.2, 4.3, 4.4


## 1.1 Systems of Linear Equations
A **linear equation** in the variables x<sub>1</sub>, ..., x<sub>n</sub> is an equation that can be written in the form a<sub>1</sub>x<sub>1</sub> + a<sub>2</sub>x<sub>2</sub> + ... + a<sub>n</sub>x<sub>n</sub> = b where *b* and the **coefficients** a<sub>1</sub>, ..., a<sub>n</sub> are real or complex numbers. The subscript *n* may be any positive integer.

A **system of linear equations** (or a **linear system**) is a collection of one or more linear equations involving the same variables - say, x<sub>1</sub>, ..., x<sub>n</sub>.

A **solution** of the system is a list (s<sub>1</sub>, s<sub>2</sub>, ..., s<sub>n</sub>) of numbers that makes each equation a true statement when the values s<sub>1</sub>, ..., s<sub>n</sub> are substituted for x<sub>1</sub>, ..., x<sub>n</sub>, respectively.

The set of all possible solutions is called the **solution set** of the linear system. Two linear systems are called **equivalent** if they have the same solution set.

A system of linear equations has
1. no solution, or
2. exactly one solution, or
3. infinitely many solutions.

A system of linear equations is said to be **consistent** if it has either one solution or infinitely many solutions; a system is **inconsistent** if it has no solution.

The essential information of a linear system can be recorded compactly in a rectangular array called a **matrix**. Given the system with the coefficients of each variable aligned in the columns, the matrix is called the **coefficient matrix** (or **matrix of coefficients**) of the system. An **augmented matrix** of a system consists of the coefficient matrix with an added column containing the constants from the right sides of the equations.

The **size** of a matrix tells how many rows and columns it has. An **m x n matrix** is a rectangular array of numbers with *m* rows and *n* columns.

**Elementary Row Operations**
1. (Replacement) Replace one row by the sum of itself and a multiple of another row.
2. (Interchange) Interchange two rows.
3. (Scaling) Multiply all entries in a row by a nonzero constant.

Two matrices are called **row equivalent** if there is a sequence of elementary row operations that transforms one matrix into the other.

If the augmented matrices of two linear systems are row equivalent, then the two systems have the same solution set.

**Two Fundamental Questions about a Linear System**
1. Is the system consistent; that is, does at least one solution *exist*?
2. If a solution exists, is it the *only* one; that is, is the solution *unique*?


## 1.2 Row Reduction and Echelon Forms

A rectangular matrix is in **echelon form** (or **row echelon form**) if it has the following three properties:
1. All nonzero rows are above any rows of all zeros.
2. Each leading entry of a row is in a column to the right of the leading entry of the row above it.
3. All entries in a column below a leading entry are zeros.

If a matrix in echelon form satisfies the following additional conditions, then it is in **reduced echelon form** (or **reduced row echelon form**):

4. The leading entry in each nonzero row is 1.
5. Each leading 1 is the only nonzero entry in its column.

Any nonzero matrix may be **row reduced** (that is, transformed by elementary row operations) into more than one matrix in echelon form, using different sequences of row operations.

**Uniqueness of the Reduced Echelon Form**\
Each matrix is row equivalent to one and only one reduced echelon matrix.

A **pivot position** in a matrix A is a location in A that corresponds to a leading 1 in the reduced echelon form of A. A **pivot column** is a column of A that contains a pivot position.

A **pivot** is a nonzero number in a pivot position that is used as needed to create zeros via row operations.

The variables corresponding to the pivot columns in the matrix are called **basic** variables. The other variables are called **free** variables.

**Existence and Uniqueness Theorem**\
A linear system is consistent if and only if the rightmost column of the augmented matrix is not a pivot column - that is, if and only if an echelon form of the augmented matrix has no row of the form [0 ... 0 b] with b nonzero. If a linear system is consistent, then the solution set contains either (i) a unique solution, when there are no free variables, or (ii) infinitely many solutions, when there is at least one free variable.


## 1.3 Vector Equations

A matrix with only one column is called a **column vector**, or simply a **vector**. Two vectors in R2 are **equal** if and only if their corresponding entries are equal.

Given two vectors *u* and *v*, their **sum** is the vector *u+v* obtained by adding corresponding entries of *u* and *v*.

Given a vector *u* and a real number *c*, the **scalar multiple** of *u* by *c* is the vector *cu* obtained by multiplying each entry in *u* by *c*.

If *n* is a positive integer, R<sup>n</sup> (read "r-n") denotes the collection of all lists of *n* real numbers, usually written as *n* x 1 column matrices.

The vector whose entries are all zero is called the **zero vector** and is denoted by 0.

**Algebraic Properties of R<sup>n</sup>**\
For all *u*, *v*, *w*, in R<sup>n</sup> and all scalars *c* and *d*:\
(i) u + v = v + u\
(ii) (u + v) + w = u + (v + w)\
(iii) u + 0 = 0 + u = u\
(iv) u + (-u) = -u + u = 0 where -u denotes (-1)u\
(v) c(u + v) = cu + cv\
(vi) (c + d)u = cu + du\
(vii) c(du) = (cd)u\
(viii) 1u = u

Given vectors v<sub>1</sub>, v<sub>2</sub>, ..., v<sub>p</sub> in R<sup>n</sup> and given scalars c<sub>1</sub>, c<sub>2</sub>, ..., c<sub>p</sub>, the vector *y* defined by y = c<sub>1</sub>v<sub>1</sub> + ... c<sub>p</sub>v<sub>p</sub> is called a **linear combination** of v<sub>1</sub>, ..., v<sub>p</sub> with weights c<sub>1</sub>, ..., c<sub>p</sub>.

A vector equation x<sub>1</sub>**a<sub>1</sub>** + x<sub>2</sub>**a<sub>2</sub>** + ... + x<sub>n</sub>**a<sub>n</sub>** = **b** has the same solution set as the linear system whose augmented matrix is [**a<sub>1</sub>** **a<sub>2</sub>** ... **a<sub>n</sub>** **b**]. In particular, **b** can be generated by a linear combination of **a<sub>1</sub>**, ..., **a<sub>n</sub>** if and only if there exists a solution to the linear system corresponding to the matrix.

If **v<sub>1</sub>**, ..., **v<sub>p</sub>** are in R<sup>n</sup>, then the set of all linear combinations of **v<sub>1</sub>**, ..., **v<sub>p</sub>** are denoted by span{**v<sub>1</sub>**, ..., **v<sub>p</sub>**} and is called the **subset of** R<sup>n</sup> **spanned** (or **generated**) **by** **v<sub>1</sub>**, ..., **v<sub>p</sub>**. That is, span{**v<sub>1</sub>**, ..., **v<sub>p</sub>**} is the collection of all vectors that can be written in the form c<sub>1</sub>**v<sub>1</sub>** + c<sub>2</sub>**v<sub>2</sub>** + ... + c<sub>p</sub>**v<sub>p</sub>** with c<sub>1</sub>, ...,  c<sub>p</sub> scalars.


## 2.3 Characterizations of Invertible Matrices

Let A and B be square matrices. If AB = I, then A and B are both **invertible**, with B = A<sup>-1</sup> and A = B<sup>01</sup>.

Let T : R<sup>n</sup> -> R<sup>n</sup> be a linear transformation and let A be the standard matrix for T. Then T is **invertible** if and only if A is an invertible matrix. In that case, the linear transformation S given by S(x) = A<sup>-1</sup>x is the unique function satisfying the equations S(T(x)) = x for all x in R<sup>n</sup> and T(S(x)) = x for all x in R<sup>n</sup>.


## 2.8 Subspaces of R<sup>n</sup>

A **subspace** of R<sup>n</sup> is any set *H* in R<sup>n</sup> that has three properties:\
a. The zero vector is in *H*.\
b. For each *u* and *v* in *H*, the sum *u+v* is in *H*.\
c. For each *u* in *H* and scalar *c*, the vector *cu* is in *H*.

The **column space** of a matrix *A* is the set Col *A* of all linear combinations of the columns of *A*.

The **null space** of a matrix *A* is the set Nul *A* of all solutions of the homogeneous equation *A*x = 0.

The null space of an *m* x *n* matrix *A* is a subspace of R<sup>n</sup>. Equivalently, the set of all solutions of a system *A*x = 0 of *m* homogeneous linear equations in *n* unknowns is a subspace of R<sup>n</sup>.

A **basis** for a subspace *H* of R<sup>n</sup> is a linearly independent set in *H* that spans *H*.

The pivot columns of a matrix *A* form a basis for the column space of *A*.


## 2.9 Dimension and Rank

Suppose the set *B* = {**b<sub>1</sub>**, ..., **b<sub>p</sub>**} is a basis for a subspace *H*. For each *x* in *H*, **the coordinates of *x* relative to the basis *B*** are weights c<sub>1</sub>, ..., c<sub>p</sub> such that x = c<sub>1</sub>b<sub>1</sub> + ... + c<sub>p</sub>b<sub>p</sub>, and the vector in R<sup>p</sup> [x]<sub>B</sub> = [c<sub>1</sub>, ..., c<sub>p</sub>] is called the **coordinate vector of *x* (relative to *B*)** or the ***B*-coordinate vector of x**.

The **dimension** of a nonzero subspace *H*, denoted by dim *H*, is the number of vectors in any basis for *H*. The dimension of the zero subspace {0} is defined to be zero.

The **rank** of a matrix *A*, denoted by rank *A*, is the dimension of the column space of *A*.

**The Rank Theorem**\
If a matrix *A* has *n* columns, then rank *A* + dim Nul *A* = *n*.

**The Basis Theorem**\
Let *H* be a *p*-dimensional subspace of R<sup>n</sup>. Any linearly independent set of exactly *p* elements in *H* is automatically a basis for *H*. Also, any set of *p* elements of *H* that spans *H* is automatically a basis for *H*.


## 3.1 Introduction to Determinants

For n >= 2, the **determinant** of an *n* x *n* matrix *A* = [a<sub>ij</sub>] is the sum of *n* terms of the form +-a<sub>1j</sub> det *A<sub>1j</sub>*, with plus and minus signs alternating, where the entries a<sub>11</sub>, a<sub>12</sub>, ..., a<sub>1n</sub> are from the first rows of *A*. In symbols, det *A* = a<sub>11</sub> det *A<sub>11</sub>* - a<sub>12</sub> det *A<sub>12</sub>* + ... + (-1)<sup>(1+n)</sup> a<sub>1n</sub> det *A<sub>1n</sub>* = sum from j = 1 to n of (-1)<sup>(1+j)</sup> a<sub>1j</sub> det *A<sub>1j</sub>*.

Given *A* = [a<sub>ij</sub>], the **(i, j) - cofactor** of A is the  number C<sub>ij</sub> given by C<sub>ij</sub> = (-1)<sup>(i+j)</sup> det *A<sub>ij</sub>*. The determinant of an *n* x *n* matrix *A* can be computed by a cofactor expansion across any row or down any column. The expansion across the *i*th row is det *A* = a<sub>i1</sub>*C<sub>i1</sub>* + a<sub>i2</sub>*C<sub>i2</sub>* + ... + a<sub>in</sub>*C<sub>in</sub>*. The expansion across the *j*th column is det *A* = a<sub>1j</sub>*C<sub>1j</sub>* + a<sub>2j</sub>*C<sub>2j</sub>* + ... + a<sub>nj</sub>*C<sub>nj</sub>*.

If *A* is a triangular matrix, then det *A* is the product of the entries on the main diagonal of *A*.


## 3.2 Properties of Determinants

**Row Operations**\
Let *A* be a square matrix.\
a. If a multiple of one row of *A* is added to another row to produce a matrix *B*, then det *B* = det *A*.\
b. If two rows of *A* are interchanged to produce *B*, then det *B* = -det *A*.\
c. If one row of *A* is multiplied by *k* to produce *B*, then det *B* = *k* x det *A*.

A square matrix *A* is invertible if and only if det A is *not* equal to 0.

If *A* is an *n* x *n* matrix, then det *A*<sup>T</sup> = det *A*.

If *A* and *B* are *n* x *n* matrices, then det(*AB*) = det(*A*) x det(*B*).


## 5.1 Eigenvectors and Eigenvalues

An **eigenvector** of an *n* x *n* matrix *A* is a nonzero vector *x* such that *Ax* = λ*x* for some scalar λ. A scalar λ is called an **eigenvalue** of *A* if there is a nontrivial solution *x* of *Ax* = λ*x*; such an *x* is called an *eigenvector corresponding to λ*.

The eigenvalues of a triangular matrix are the entries on its main diagonal.

If **v<sub>1</sub>**, ..., **v<sub>r</sub>** are eigenvectors that correspond to distinct eigenvalues λ<sub>1</sub>, ..., λ<sub>r</sub> of an *n* x *n* matrix *A*, then the set {**v<sub>1</sub>**, ..., **v<sub>r</sub>**} is linearly independent.


## 5.2 The Characteristic Equation

Let *A* be an *n* x *n* matrix, let *U* be any echelon form obtained from *A* by row replacements and row interchanges (without scaling), and let *r* be the number of such row interchanges. Then the **determinant** of *A*, written as det *A*, is (-1)<sup>r</sup> times the product of the diagonal entries *u<sub>11</sub>*, ..., *u<sub>nn</sub>*. If A is invertible, then the determinant is nonzero because the diagonal entries are all pivots.

**Properties of Determinants**\
Let *A* and *B* be *n* x *n* matrices.\
a. *A* is invertible if and only if det *A* is not 0.\
b. det *AB* = (det *A*)(det *B*).\
c. det *A<sup>T</sup>* = det *A*.\
d. If *A* is triangular, then det *A* is the product of the entries on the main diagonal of *A*.\
e. A row replacement operation on *A* does not change the determinant. A row interchange changes the sign of the determinant. A row scaling also scales the determinant by the same scalar factor.

A scalar λ is an eigenvalue of an *n* x *n* matrix *A* if and only if λ satisfies the characteristic equation det(*A*-λI) = 0.

If *A* is an *n* x *n* matrix, then det(*A*-λI) is a polynomial of degree *n* called the **characteristic polynomial** of *A*. In general, the **multiplicity** of an eigenvalue λ is its multiplicity as a root of the characteristic equation.

If *A* and *B* are *n* x *n* matrices, then *A* is **similar** to *B* if there is an invertible matrix *P* such that *P<sup>-1</sup>AP = B*, or, equivalently, *A = PBP<sup>-1</sup>*. Changing *A* into *P<sup>-1</sup>AP* is called a **similarity transformation**.

If *n* x *n* matrices *A* and *B* are similar, then they have the same characteristic polynomial and hence the same eigenvalues(with the same multiplicities).


## 5.3 Diagonalization

A square matrix *A* is said to be **diagonalizable** if *A* is similar to a diagonal matrix, that is, if *A = PDP<sup>-1</sup>* for some invertible matrix *P* and some diagonal matrix *D*.

**The Diagonalization Theorem**\
An *n* x *n* matrix *A* is diagonalizable if and only if *A* has *n* linearly independent eigenvectors. In fact, *A = PDP<sup>-1</sup>*, with *D* a diagonal matrix, if and only if the columns of *P* are *n* linearly independent eigenvectors of *A*. In this case, the diagonal entries of *D* are eigenvalues of *A* the correspond to, respectively, to the eigenvectors in *P*. In otherwords, *A* is diagonalizable if and only if there are enough eigenvectors to form a basis for R<sup>n</sup>. We call such a basis an **eigenvector basis** of R<sup>n</sup>.

An *n* x *n* matrix with *n* distinct eigenvalues is diagonalizable.

Let *A* be an *n* x *n* matrix whose distinct eigenvalues are λ<sub>1</sub>, ..., λ<sub>p</sub>.\
a. For 1 <= *k* <= *p*, the dimensions of the eigenspace for λ<sub>k</sub> is less than or equal to the multiplicity of he eigenvalue λ<sub>k</sub>.\
b. The matrix *A* is diagonalizable if and only if the sum of the dimensions of the eigenspaces equals *n*, and this happens if and only if (i) the characteristic polynomial factors completely into linear factors and (ii) the dimension of the eigenspace for each λ<sub>k</sub> equals the multiplicity of λ<sub>k</sub>.\
c. If *A* is diagonalizable and *B<sub>k</sub>* is a basis for the eigenspace corresponding to λ<sub>k</sub> for each *k*, then the total collection of vectors in the set *B<sub>1</sub>*, ..., *B<sub>p</sub>* forms an eigenvector basis for R<sup>n</sup>.


## 4.1 Vector Spaces and Subspaces

A **vector space** is a nonempty set *V* of objects, called *vectors*, on which are defined two operations, called *addition* and *multiplication by scalars*, subject to 10 axioms (or rules) listed below. The axioms must hold for all vectors **u**, **v**, and **w** in *V* and for all scalars *c* and *d*.
1. The sum of *u* and *v*, denoted by *u+v*, is in *V*.
2. *u* + *v* = *v* + *u*.
3. (*u* + *v*) + *w* = *u* + (*v* + *w*).
4. There is a zero vector 0 in *V* such that *u* + 0 = *u*.
5. For each *u* in *V*, there is a vector *-u* in *V* such that *u* + (*-u*) = 0.
6. The scalar multiple of *u* by *c*, denoted by *cu*, is in *V*.
7. *c*(*u* + *v*) = *cu* + *cv*.
8. (*c* + *d*)*u* = *cu* + *du*.
9. *c*(*du*) = (*cd*)*u*.
10. 1*u* = *u*.

A **subspace** of a vector space *V* is a subset *H* of *V* that has three properties:\
a. The zero vector of *V* is in *H*.\
b. *H* is closed under vector addition.\
c. *H* is closed under scalar multiplication.

If **v<sub>1</sub>**, ..., **v<sub>p</sub>** are in a vector space *V*, then span{**v<sub>1</sub>**, ..., **v<sub>p</sub>**} is a subspace of *V*.

## 4.2 Null Spaces, Column Spaces, and Linear Transformations

The **null space** of an *m* x *n* matrix *A*, written as Nul *A*, is the set of all solutions of the homogeneous equation *Ax* = 0. In set notation, Nul *A* = {*x*: *x* in R<sup>n</sup> and *Ax* = 0}. The null space is a subspace of R<sup>n</sup>. Equivalently, the set of all solutions of a system *A*x = 0 of *m* homogeneous linear equations in *n* unknowns is a subspace of R<sup>n</sup>.

The **column space** of an *m* x *n* matrix *A*, written as Col *A*, is the set of all linear combinations of the columns of *A*. If *A* = [**a<sub>1</sub>** ... **a<sub>n</sub>**], then Col *A* = span{**a<sub>1</sub>**, ..., **a<sub>n</sub>**}. The column spaces is a subspace of R<sup>m</sup>. The column space is all of R<sup>m</sup> if and only if the equation *Ax* = *b* has a solution for each *b* in R<sup>m</sup>.

A **linear transformation** *T* from a vector space *V* to a vector space *W* is a rule that assigns to each vector *x* in *V* a unique vector *T(x)* in *W*, such that (i) *T(u + v)* = *T(u)* + *T(v)* for all *u*, *v*, in *V*, and (ii) *T(cu)* = *cT(u)* for all *u* in *V* and all scalars *c*.

The **kernal** (or **null space**) of such a *T* is the set of all *u* in *V* such that *T(u)* = 0.

The **range** of *T* is the set of all vectors in *W* of the form *T(x)* for some *x* in *V*.


## 4.3 Linearly Independent Sets; Bases

An indexed set {**v<sub>1</sub>**, ..., **v<sub>p</sub>**} of two or more vectors, with **v<sub>1</sub>** not equal to 0 is lnearly dependent if and only if some **v<sub>j</sub>** (with *j* > 1) is a linear combination of the preceding vectors, **v<sub>1</sub>**, ..., **v<sub>j-1</sub>**.

Let *H* be a subspace of vector space *V*. An indexed set of vectors *B* = {**b<sub>1</sub>**, ..., **b<sub>p</sub>**} in *V* is a **basis** for *H* if (i) *B* is a linearly independent set, and (ii) the subspace spanned by *B* coincides with *H*; that is, *H* = span{**b<sub>1</sub>**, ..., **b<sub>p</sub>**}.

**The Spanning Set Theorem**\
Let *S* = {**v<sub>1</sub>**, ..., **v<sub>p</sub>**} be a set in *V*, and let *H* = span{**v<sub>1</sub>**, ..., **v<sub>p<sub>**}.\
a. If one of the vectors in *S* - say, **v<sub>k</sub>** - is a linear combination of the remaining vectors in *S*, then the set formed from *S* by removing **v<sub>k</sub>** still spans *H*.\
b. IF *H* not equal to {0}, some subset of *S* is a basis for *H*.

The pivot columns of a matrix *A* form a basis for Col *A*.


## 4.4 Coordinate Systems

**The Unique Representation Theorem**\
Let *B* = {**b<sub>1</sub>**, ..., **b<sub>n</sub>**} be a basis for a vector space *V*. Then for each *x* in *V*, there exists a unique set of scalars c<sub>1</sub>, ..., c<sub>n</sub> such that *x* = c<sub>1</sub>**b<sub>1</sub>** + ... + c<sub>n</sub>**b<sub>n</sub>**. The **coordinates of x relative to the basis *B*** (or the ***B*-coordinates of x**) are the weights c<sub>1</sub>, ..., c<sub>n</sub>.

Let *B* = {**b<sub>1</sub>**, ..., **b<sub>n</sub>**} be a basis for a vector space *V*. Then the coordinate mapping *x* -> [*x*]<sub>*B*</sub> is a one-to-one linear transformation from *V* onto R<sub>n</sub>.

In general, a one-to-one linear transformation from a vector space *V* to a vector space *W* is called an **isomorphism** from *V* to *W*.

## Overarching

**The Invertible Matrix Theorem**\
Let *A* be a *n* x *n* matrix. Then, the following statements are equivalent.

a. *A* is an invertible matrix.\
b. *A* is row equivalent to the *n* x *n* identity matrix.\
c. *A* has *n* pivot positions.\
d. The equation *A*x = 0 has only the trivial solution.\
e. The columns of *A* form a linearly independent set.\
f. The linear transformation x -> *A*x is one-to-one.\
g. The equation *A*x = b has at least one solution for each b in R<sup>n</sup>.\
h. The columns of *A* span R<sup>n</sup>.\
i. The linear transformation x -> *A*x maps R<sup>n</sup> onto R<sup>n</sup>.\
j. There is an *n* x *n* matrix C such that C*A* = I.\
k. There is an *n* x *n* matrix D such that *A*D = I.\
l. *A*<sup>T</sup> is an invertible matrix.\
m. The columns of *A* form a basis of R<sup>n</sup>.\
n. Col *A* = R<sup>n</sup>.\
o. dim Col *A* = n.\
p. rank *A* = n.\
q. Nul *A* = {0}.\
r. dim Nul *A* = 0.\
s. The number 0 is *not* an eigenvalue of *A*.\
t. The determinant of *A* is *not* 0.
