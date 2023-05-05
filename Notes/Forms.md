# Forms

Some definitions about the fundamental objects that will help us to, among others, calculate distances and fully define integration in curved spaces.

----

## Covectors

**<u>Definition:</u>** Given a vector space $V$ over a field $\mathbb{K} \in \{\mathbb{R},\mathbb{C}\}$, the **dual space** $V^\wedge$ is the space of all linear maps $a:V\to \mathbb{K}$ , i.e. for any $u,v \in V,$ and $\lambda \in \mathbb{K}$
$$
a(v + \lambda u) = a(v) + \lambda a(u) 
$$
**<u>*Proposition:*</u>** The dual space $V^\wedge$ is isomorphic to $V$. 

**Proof:** Consider a basis $\{v_i\}$ on $V$ we will show that we can find a basis on $V^\wedge$ related to $\{v_i\}$ . Consider an element $a \in V^\wedge$ and $u \in V$ where $u = u^i v_i$.  By definition we know that
$$
a(u) = a\left(\sum_{i=1}^{\dim V}u^iv_i\right) = \sum_{i=1}^{\dim V}u^i a(v_i)
$$
Take the vectors $\{v^1,v^2,\cdots,v^n\} \sub V^\wedge$ to be defined as follows
$$
v^i(v_j) = \delta^i_j
$$
$v^i(u)=v^i(u^jv_j)=u^jv^i(v_j)=u^j\delta^i_j=u^i$

Take any $a \in V^\wedge$ 
$$
a(u)=u^ia(v_i)=v^i(u)a(v_i) = a_iv^i(u)
$$
Since this is true for any $u$ we can write $a \in V^\wedge$ as
$$
a=a_iv^i
$$
Proving that our set of covectors <u>spans</u> $V^\wedge$. To show that they are linearly indepentent we need to show that
$$
a=a_iv^i=0 \in V^\wedge \iff  a_i = 0 \in \mathbb{K}
$$
Forward if $a=0$ we have $0 = a(v_i)= a_jv^j(v_i)  = a_i$. The backward is trivial if $a_i = 0$ then $a(u) = v^i(u)\cdot0 = 0$  $\Box$

**<u>*Corollary:*</u>** Just like we can think of dual vectors as maps of vectors, we can think of vectors as maps of dual vectors, i.e. $(V^\wedge)^\wedge \cong V$.



## Covectors on Manifolds

We have already discussed the construction of a vector space at every point in a manifold. This natural vector space at a point $p$ on a manifold $M$ is called the tangent space $T_p M$ and it contains vectors that we additionally interpret as derivations at that point. Let’s see what happens at the covector space.

**<u>Definition:</u>** Given a smooth manifold $M$ and a point $p \in M$ we define the dual space of the tangent space $T_pM$ as
$$
T_p^*M = (T_pM)^\wedge
$$
called the **cotangent space**. We additionally define the **cotangent bundle** as the disjoint union of the cotangent spaces like so:
$$
T^*M=\coprod_{p\in M} T^*_pM
$$


To uncover the differential nature of these objects, consider the definition of the differential below

**<u>Definition:</u>** Given a function $f \in C^\infty(M)$ we can define its **differential** at $p\in M$ as a linear map of vectors defined like so
$$
df_p(X_p) = X_pf,\ \forall X_p  \in T_pM
$$
**<u>Corollary:</u>** The differential $df_p \in T_p^*M$ is an element of the cotangent space at that point, as a linear map of vectors. 



**<u>*Proposition:*</u>** Given some point $p\in M$ and some coordinates $\{x^i\}$ in a neihborhood around that point we can define a basis for the cotangent space by taking their differentials $\{dx^i_p\}$.

*Proof:* Consider the standard basis of the tangent space $T_pM$ $\{\left.\frac{\partial }{\partial x^i} \right|_p\}$. We can show that:
$$
dx^i\left(\left.\frac{\partial }{\partial x^j} \right|_p \right) = \left.\frac{\partial x^i}{\partial x^j} \right|_p= \delta^i_j
$$
Now consider an arbitrary covector $\alpha \in T^*_pM$ and a vector $X\in T_pM$.
$$
\alpha(X) = \alpha\left(X^i \left.\frac{\part }{\part x^i}\right|\right) = X^i \alpha\left(\left.\frac{\part }{\part x^i}\right|\right) = dx^i(X)\  \alpha\left(\left.\frac{\part }{\part x^i}\right|\right) = a_i\, dx^i(X)
$$
The last step was carried out using the property we proved in (10). This proves that we can write every covector in this basis. With the same equation we have proved bijectivity so it is a dual basis. $\Box$



Similarly to vector fields we can have covector fields where they assign a covector at every point in the manifold. We call those 1-forms

**<u>Definition:</u>** A covector field or a **1-form** $\omega$ is a map that assigns at each point $p \in M$ an element of the cotangent space at that point
$$
\begin{align*}
\omega : M &\to T^*M\\
p &\mapsto \omega(p) = \omega_p \in T_p^*M \subset T^*M
\end{align*}
$$
We denote the space of all one forms on manifold $M$ by $\Omega^1(M)$.

**<u>Corollary:</u>** Notice that the 1-form $df$ is the covector fields that evaluates to $df(p) = df_p$! Those differentials are the one forms!



## Tensors

One forms and vectors are tensors. Here we will quickly define tensors, what is the order of a tensor, and a bunch of products that can help us construct higher order tensors.



**<u>Definition:</u>** A **tensor** $T$ **of** **order** $k$ is a $k-$linear map from $k$ copies of the vector space $V$ or its dual $v^\vee$ to the underlying field (in this case the real numbers). In particular if we have two positive integers $n,m$ such that $n+m=k$ then a $(n,m)$ tensor $T$ is defined as a $k$-linear map
$$
T:V^n\times \left(V^\vee \right)^m \to \mathbb{R}
$$
We define the set of all $(n,m)$ tensors as $\mathcal{T}^{n,m}$.



We can clearly see that vectors in this notation are $(0,1)$ tensors, while covectors are $(1,0)$ tensors. But we can have more! For example the determinant of $\mathbb{R}^n$ $\det$ is a $(n,0)$ tensor! 

*<u>Example:</u>* Thedeterminant $\det$ of $\mathbb{R}^n$ is a $(n,0)$ tensor.

Notice that for $\mathbb{R}^n$ the determinant is a map that assigns a real number to an $n \times n$ matrix. But such a matrix is just a collection of $n$ (column) vectors. In that notation we can write the determinant as:
$$
\det(v_1,v_2,\cdots,v_n) = \varepsilon_{i_1, i_2,\cdots, i_n}v_1^{i_1} v_2^{i_2} \cdots v_n^{i_n}
$$
With this definition this is clearly a $k$-linear map of vectors and this an $(n,0)$ alternating tensor. 



But how can I create new tensors? You ask. Lemme tell you! 

**<u>Defintion:</u>** The **tensor product** of a $k$ tensor $a$ and an $l$ tensor $b$ is defined as the $k+l$ tensor $a\otimes b$ that acts like so:
$$
a\otimes b\ (v_1,v_2,\cdots,v_{k+l}) = a(v_1,v_2,\cdots,v_k)\ b(v_{k+1},v_{k+2},\cdots,v_{k+l})
$$
This is an eazy way to create higher order tensors from lower order ones, but we can use it to also create other types of products taht produce tensors with certain properties. 



## K-Forms

We have already seen what a one-form is as a covector field. Here we will explicitly talk about $k$-forms on a manifold.

**<u>Definition:</u>** Given a smooth manifold $M$, a $k$-form $\omega$ is an antisymmetric $(k,0)$ tensor field over the manifold $M$ (that is that the $k$-tensor takes $k$ vectors from the tangent bundle of M). We call the space of $k$-forms $\Omega^k(M)$. To be more precise about the antisymmetricity of $\omega$ consider a permutation $\sigma$, then

$$
\omega(v_1,v_2,\cdots,v_n) = \text{sign}(\sigma)\ \omega(v_{\sigma(1)},v_{\sigma(2)},\cdots,v_{\sigma(n)})
$$


The next thing we want to define is a product of k-forms that keeps their antisymmetricity. This product is called the wedge product

**<u>Definition:</u>** Given $\omega \in \Omega^k(M)$ and $\eta \in \Omega^l(M)$, their wedge product is given by
$$
\begin{align*}
\cdot \wedge \cdot:\Omega^k(M)\times \Omega^l(M) &\to \Omega^{k+l}(M)\\
(\omega,\eta) &\mapsto\omega\wedge \eta
\end{align*}
$$
This is given precisely by this formula
$$
\omega \wedge \eta\ (v_1,v_2,\cdots,v_{k+l}) = \sum_{\sigma \in S_{k+l}}\frac{1}{k!\,l!} \omega(v_{\sigma(1)},v_{\sigma(2)},\cdots,v_{\sigma(k)})\ \eta(v_{\sigma(k+1)},v_{\sigma(k+2)},\cdots,v_{\sigma(k+l)})
$$


Ok ok ok we are getting carried away in definitions that seem too disconnected from what we were talking about previously! We said 1-forms are the $df$ stuff! What is going on with this new way of thikning about them? Well clearly they must be compatible, but using the basis of 1 forms we can create a basis for the k-forms



**<u>Proposition:</u>** Given $\{dx^i\}$ the standard basis of $\Omega^1(M)$, the wedge products of the basis form a basis for the $k$-forms $\{dx^{i_1} \wedge dx^{i_2} \wedge \cdots \wedge dx^{i_k} \}$.

The proof is by calculation, but it’s just tedius, so I won’t do it. But this leads to a very interesting corollary.

**<u>Corollary:</u>** On an $n$ dimensional manifold $M$ we can define $k$-forms with $0\leq k\leq n$. 

In other words $\Omega^{n+1}(M) = \empty$.

*Proof:* We will show that considering $\Omega^k(M)$ as vector spaces over $C^\infty(M)$ their dimensions follow pascall’s triangle. 

First consider a direct property of the antisymetricity of forms. Given any form $\omega \in \Omega^k(M)$ we have that
$$
\omega\wedge \omega = - \omega\wedge\omega = 0
$$
Now take a look at the basis for $\Omega^k(M)$.
$$
dx^{i_1}\wedge dx^{i_2} \wedge \cdots \wedge dx^{i_k}
$$
if for some $i_a = i_b$ then this vector is $0$ so it does not belong to the basis. For example $dx^\wedge dx = 0$ and it does not belong to the basis for $\Omega^2(M)$. Therefore, each basis vector is composed of the wedge product of $k$ different $dx^i$s. The number of basis vectors of $\Omega^k(M)$ is the number ways we can pick $k$ elements out of $n$ without replacement. In other word
$$
\dim\Omega^k(M) = \begin{pmatrix}n\\k\end{pmatrix}
$$
Notice that if we try to construct a form out of the wedge product of more than $n$ basis, we are guaranteed to have repeated entries, so we cant have forms with $k>n$. $\Box$



## The exterior derivative

The exterior derivative is a way to create higher order forms. We have been silently using it, but here is a formal definition. 



**<u>Definition:</u>** Given a smooth manifold $M$, the **exterior derivative** $d$ is a map
$$
d : \Omega^k(M) \to \Omega^{k+1} (M)
$$
That is defined like so. For $f \in C^\infty(M)$
$$
df = \frac{\partial f}{\part x^i} dx^i
$$
(Note that this is compatible with our previous, pointwise, definition of the differential)

For $\omega \in \Omega^k(M)$ this is defined like
$$
\begin{align*}
d\omega 
&= d\omega_i \wedge dx^i\\
&= \frac{\part \omega_i}{\part x^j}dx^j \wedge dx^i

\end{align*}
$$
**<u>Proposition:</u>** Given two forms $\omega \in \Omega^k(M)$ and $y \in \Omega^l(M)$ then $d(\omega \wedge \eta ) = d\omega \wedge\eta + \omega \wedge d\eta$.

This is proved by applying the definition. BTW this is product rule!



***Note:*** Check this out! This means that we can identify $C^\infty(M) = \Omega^0(M)$ as the zero-forms on $M$, since $df \in \Omega^1(M)\ \forall f\in C^\infty(M)$. So that’s nice, it completes the pascal triangle. 











