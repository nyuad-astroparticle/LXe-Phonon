# Automorphism Groups of Vector Spaces

Vector spaces are really cool, and thinking of formalizing the way we can linearly move vectors from one place to another gives rize to many beautiful Lie Groups. Here is a survey of them because practice makes perfect.

-----

## Graded Alebras

Algebras are vectorspaces with extra structure that can almost resemble a field, but without the requirement of explicit inverses. There are 3 infinite graded division algebras by *Cartan’s Theorem.* These are: $\mathbb{R}, \mathbb{C}, \text{ and } \mathbb{H}$. We can create vector spaces over them and imbude them with the following inner products we are going to call **canonical** because of their relation to the eucledian distance.

**<u>Definition:</u>** The following are the **canonical inner products** defined on vector spaces over the corresponding graded division algebras

1. In the **real numbers** $\mathbb{R}$ we have the **Eucledian Scalar Product**:
   $$
   \lang u,v \rang \coloneqq u^Tv = u_i\, v^i
   $$

2. In the **complex numbers** $\mathbb{C}$ we have the **Hermitian Scalar Product**:
   $$
   \lang u , v \rang \coloneqq u^\dagger v = \bar{u}_i\, v^i
   $$

3. In the **quarternions** $\mathbb{H}$ we have the **Symplectic Scalar Product**:
   $$
   \lang  u, v \rang \coloneqq u^\dagger v = \bar{u}_i v^i
   $$

**<u>Definition:</u>** The canonical volume form on vector space $V$ over an algebra $\mathbb{K} \in \{\mathbb{R},\mathbb{C},\mathbb{H}\}$ is known as the **determinant** and is defined by
$$
\begin{align}
\text{vol}: V\times \cdots\times V &\to \mathbb{K}\\
(v_1,v_2,\cdots,v_n) &\mapsto \det(v_1,v_2,\cdots,v_n)
\end{align}
$$
where $v_i \in V$ is a vector. 

## The Linear Automorphism Groups

We have our vectorspaces and their canonical notion of distance (using their inner products). What we need to do now is to classify their automorphisms. As it turns out they can form some extremely useful groups.



**<u>Definition:</u>** Here are some Linear Automorphism Groups (under composition) of vector spaces over $\mathbb{K} \in \{\mathbb{R},\mathbb{C},\mathbb{H}\}$

1. The **general linear group** of order $n$ is defined as follows
   $$
   \text{GL}(n,\mathbb{K}) = \{A: \mathbb{K}^n \to \mathbb{K}^n\ \mid \ A \text{ linear and invertible}\}
   $$
   This can further be specialized as $A$ being an $n\times n$ matrix with entries in $\mathbb{K}$ and $\det{A} \neq 0$.

2. The next one of interest is the **special linear group** of order $n$ defined like so
   $$
   \text{SL}(n,\mathbb{K})=\{A \in \text{GL}(n,\mathbb{K})\, \mid \, \text{vol}(Av_1,Av_2,\cdots,Av_n) = \text{vol}(v_1,v_2,\cdots,v_n)\}
   $$
   This requirement can be boiled down to the statement that $\det A=1$.



After these two big groups we can further define groups that originate from norm preserving linear automorphisms like so

1. The **Orthogonal Group** of order $n$ is given by
   $$
   \text{O}(n) = \{A\in \text{GL}(n,\mathbb{R})\, \mid \, \lang Au,Av \rang = \lang u,v\rang \ \forall u,v \in \mathbb{R}^n\}
   $$

2. The **Unitary Group** of order $n$ is given by
   $$
   \text{U}(n) = \{A\in \text{GL}(n,\mathbb{C})\, \mid \, \lang Au,Av \rang = \lang u,v\rang \ \forall u,v \in \mathbb{C}^n\}
   $$

3. The **Symplectic Group** of order $n$ is given by
   $$
   \text{Sp}(n) = \{A\in \text{GL}(n,\mathbb{H})\, \mid \, \lang Au,Av \rang = \lang u,v\rang \ \forall u,v \in \mathbb{H}^n\}
   $$

4. 



Similarly each once can get the prefix **Special** by simply intersecting it with the special linear group of order $n$ in order to consider only linear automorphisms that preserve the volume form.

1. The **Special Orthogonal Group** of order $n$ 
   $$
   \text{SO}(n)= \text{O}(n) \cap \text{SL}(n,\mathbb{R})
   $$

2. The **Special Unitary Group** of order $n$ 
   $$
   \text{SU}(n)= \text{U}(n) \cap \text{SL}(n,\mathbb{C})
   $$

These sets, under composition, are called the **classical groups** for some reason.



**<u>*Proposition:*</u>** (Alternative descriptions of Classical Groups) For $n\geq 1$ we have

1. The special linear groups are given by
   $$
   \text{SL}(n,\mathbb{K}) = \{A \in \text{GL}(n,\mathbb{K})\,\mid\, \det(A) = 1\}
   $$

2. The orthogonal, unitary, and symplectic groups can be given by
   $$
   \begin{align}
   \text{O}(n) &= \{A \in \text{GL}(n,\mathbb{R})\,\mid\, A A^T = I\}\\
   \text{U}(n) &= \{A \in \text{GL}(n,\mathbb{C})\,\mid\, A A^\dagger = I\}\\
   \text{O}(n) &= \{A \in \text{GL}(n,\mathbb{H})\,\mid\, A A^\dagger = I\}\\
   \end{align}
   $$

3. The special orthogonal and special unitary groups can be given by
   $$
   \begin{align}
   \text{SO}(n) &= \{A \in \text{GL}(n,\mathbb{R})\,\mid\, A A^T = I,\ \det(A) = 1\}\\
   \text{SU}(n) &= \{A \in \text{GL}(n,\mathbb{C})\,\mid\, A A^\dagger = I,\ \det(A) = 1\}\\
   \end{align}
   $$



<u>*Note:*</u> Elements of the symplectic group are already volume preserving i.e. $\text{Sp}(n) \subset \text{SL}(n,\mathbb{H})$.



**<u>Theorem:</u>** As subgroups of general linear groups, classical groups are linear. They also have the following dimensions as vector spaces
$$
\begin{align*}
\dim \text{SL}(n,\mathbb{R}) &= n^2 - 1\\
\dim \text{SL}(n,\mathbb{C}) &= 2n^2 - 2\\
\dim \text{SL}(n,\mathbb{H}) &= 4n^2 - 1\\
\\
\dim \text{O}(n) &= \frac{1}{2}n(n-1)\\
\dim \text{U}(n) &= n^2\\
\dim \text{Sp}(n) &= 2n^2 + n\\
\\
\dim \text{SO}(n) &= \frac{1}{2}n(n-1)\\
\dim \text{SU}(n) &= n^2 - 1\\
\end{align*}
$$


## Lie Algebras of Classical Groups

Along with knowing the groups, it is nice to know their lie Algebras. Before we do that take a look at this theorem! It’s so cool!

**<u>Theorem: (Lie’s Third Theorem)</u>** Every finite dimensional real Lie algebra is isomorphic to the Lie algebra of some connected Lie group.

This just says that if I give you a Lie algebra, I can make a shape that has this algebra. In other words imposing the symmetries in your space using a Lie algebra, defines its shape.

**<u>Theorem:</u>** These are the Lie algebras of the Classical Groups we discussed above
$$
\begin{align*}
\mathfrak{gl}(n,\mathbb{K}) &= \text{Mat}(n\times n, \mathbb{K})\\
\\
\mathfrak{sl}(n,\mathbb{R}) &= \{A \in \text{Mat}(n\times n,\mathbb{R})\,\mid\, \tr(A)=0\}\\

\mathfrak{sl}(n,\mathbb{c}) &= \{A \in \text{Mat}(n\times n,\mathbb{C})\,\mid\, \tr(A)=0\}\\

\mathfrak{sl}(n,\mathbb{H}) &= \{A \in \text{Mat}(n\times n,\mathbb{H})\,\mid\, \Re(\tr(A))=0\}\\
\\
\mathfrak{o}(n) &= \{A \in \text{Mat}(n\times n,\mathbb{R})\,\mid\, A + A^T = 0\}\\

\mathfrak{u}(n) &= \{A \in \text{Mat}(n\times n,\mathbb{C})\,\mid\, A + A^\dagger = 0\}\\

\mathfrak{sp}(n) &= \{A \in \text{Mat}(n\times n,\mathbb{H})\,\mid\, A + A^\dagger = 0\}\\
\\
\mathfrak{so}(n) &= \mathfrak{o}(n)\\

\mathfrak{su}(n) &= \{A \in \text{Mat}(n\times n,\mathbb{C})\,\mid\, A + A^\dagger = 0,\ \tr(A)=0\}
\end{align*}
$$

## Examples

Let’s do a physics example just so that we can actually be doing something meaningful. The **Pauli matrices**. The Lie algebra $\mathfrak{su}(2)$ has dimension 3 and consists of the $2\times 2$ skew-hermitian matrices of trace $0$. 

To find the elements here are three linearly independent elements
$$
\begin{align*}
\sigma_1 = \begin{pmatrix} 0 & 1 \\ 1 & 0\end{pmatrix} && \sigma_2 = \begin{pmatrix} 0 & -i \\ i & 0\end{pmatrix} && \sigma_3 = \begin{pmatrix} 1 & 0 \\ 0 & -1\end{pmatrix}
\end{align*}
$$
These definitely form a basis for $\mathfrak{su}(2)$ but we can get a better one defined like so
$$
\tau_a = -\frac{i}{2}\sigma_a
$$
Why is this nice? Because it has a pretty commutation relation
$$
[\tau_a,\tau_b]=\epsilon_{abc}\tau_c
$$
We can map now the bases of $\mathfrak{so}(3)$ to the bases of $\mathfrak{su}(2)$ and discover that we can create a Lie algebra isomorphism!



















