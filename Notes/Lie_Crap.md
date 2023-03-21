# Lie Crap

Some notes on the definitions of structures that are useful in describing groups that are also manifolds.

## Lie Groups & Algebras

The best place to start is by defining a Lie Group and its corresponding Lie algebra.

**<u>Definition:</u>** A Lie Group G is a group in the Algebraic Sense, and a smooth manifold with a multiplication map
$$
\mu:G\times G\to G\\
(a,b)\mapsto a\cdot b
$$
And an inversion map
$$
i:G\to G\\
a\mapsto a^{-1}
$$
that are both smooth as maps between manifolds.



**<u>Definition:</u>** A real Algebra A is a finite dimensional real vector space $(A,\mathbb{R})$ with a bilinear map
$$
\cdot:A\times A \to A\\
(a,b) \mapsto a\cdot b
$$
and a unit element $1\in A$ such that $1\cdot a = a = a \cdot 1,\ \forall a \in A$.

- $A$ is <u>normed</u> iff $\exists $ a norm $||\cdot||$ on $A$ such that:
  $$
  ||\cdot|| : A \to A\\
  ||a\cdot b|| = ||a|| \cdot ||b||,\ \forall a,b \in A
  $$

- $A$ is a <u>division</u> algebra iff the following is true:
  $$
  a\cdot b = 0 \implies a=0 \text{ or } b = 0
  $$

**<u>Definition:</u>** A Lie Algebra $A$ is a vector space $\mathbb{V}$ with an inner product called the **Lie Bracket** given by:
$$
[\cdot,\cdot]:\mathbb{V}\times \mathbb{V}\to \mathbb{V}
$$
such that:

1. $[\cdot,\cdot]$ is bilinear
2. $[a,b] = - [b,a]\ \forall a,b\in \mathbb{V}$
3. $[a,[b,c]] + [c,[a,b]] + [b,[c,a]] = 0,\ \forall a,b,c\in \mathbb{V}$ (Jacobi Identity)



Lie groups can be used to form multiple Lie Algebras, however, there is a particular Lie algebra that is very useful. So useful in fact that we associate it with the group by calling it *it’s Lie Algebra.* 

**<u>Definition:</u>** Given a Lie Group G, and two vector fields $X,Y \in \mathfrak{X}(G)$ we define the **commutator** of vector fields $[\cdot,\cdot]$ like so
$$
[\cdot,\cdot]:\mathfrak{X}(G)\times\mathfrak{X}(G)\to\mathfrak{X}(G)\\
(X,Y)\mapsto[X,Y] \coloneqq X\circ Y - Y \circ X
$$
**<u>*Corollary:*</u>** The space of vector fields on $G$, $\mathfrak{X}(G)$, together with the commutator of vector fields, is a Lie Algebra.

---

Small diversion on pushforward of vector fields. The pushforward is normally define for vectors form one tangent space to another given a smooth map of manifolds. However, that map being smooth does not guarantee that pushing forward each smooth vector of a vector field would produce a smooth vector field on the target manifold. Let’s formalize this using a *related vector fields*.

**<u>Definition:</u>** Given a smooth map of manifolds $F:M\to N$, two vector fields $X\in \mathfrak{X}(M)$ and $Y\in \mathfrak{X}(N)$ are related iff
$$
F_{*,p}X=Y_{F(p)},\ \forall p\in M
$$
Notice that with this definition we don’t run into invertibility issues. Also notice that there is no restriction that the vector field $Y$ is unique. This would happen only iff $F$ is a diffeomorphism of manifolds. 



Let’s now formalize some stupid notation that can take one ages to fucking realize…

**<u>Notation:</u>** Given a smooth map of manifolds $F:M\to N$ we can find it’s pushforward as the map
$$
F_*:TM\to TN\\
x\in TM \to F_*x\in TN = F_{*,\pi_M(x)}x
$$
Furthermore, if $F$ is a diffeomorphism we can push vector fields of $M$ as sections of $TM$. Given $X\in \mathfrak{X}(M)$ we have
$$
F_*X=Y\in \mathfrak{X}(N)\\
Y(p) = F_*X(p)=F_*X_{F^{-1}p}=F_{*,F^{-1}p}X_{F^{-1}p}\ \forall p\in N
$$
This notation confused me for roughly 8 hours straight!

Anyway, moving on………..

---

**<u>Definition:</u>** A vector field $X\in \mathfrak{X}(G)$ is called **left invariant** iff it is invariant under the pushforward by the left translation map $L_g$ given for some $g\in G$ by
$$
L_g:G\to G\\
a\mapsto g\cdot a.
$$
Speficially we have that for all test functions $f\in C^\infty(G)$ 
$$
L_{g*}(X) = X
$$
Notice that we can do that because $X$ is invertible

**<u>*Proposition:*</u>** The set of all left invariant vector fields on $G$ together with the commutator of vector fields, fors a Lie Algebra.

**Proof:** We know that vector fields over $G$ under addition of vector fields form a vector space. We first show that the left invariant vector fields are a subspace of $\mathfrak{X}(G)$. 

Given two vector fields $X,Y$ that are left invariant, $X+Y$ is left invariant since 
$$
L_{g*}(X+Y) = L_{g*}(X)+L_{g*}(Y) = X+Y
$$
Clearly, $0\in \mathfrak{X}(G)$ is left invariant. So the left invariant vector fields form a subspace of $\mathfrak{X}$(G) $\Box$

**<u>Definition:</u>** The set of all left invariant vector fields denoted by $\mathfrak{g}$ is known as the **Lie algebra of**  (or associated to) $G$.

**<u>*Proposition:*</u>** The Lie algebra $\mathfrak{g}$ of a Lie group $G$ is isomorphic to the tangent space at the identity $e\in G$, $T_eG$.
$$
\mathfrak{g} \cong T_eG
$$
**Proof:** The proof is based on the fact that we only need one vector to fully define a left invariant vector field on $G$. We first show that $L_{g,*} X_h=X_{gh}$ for any $g,h \in G$. This is true, we just need to apply notation, i.e.
$$
\begin{align}
(L_{g,*}X)(gh) 
&=L_{g,*} X_{L_g^{-1}(h)}\\
&=L_{g,*} X_{h}\\
&=X(gh)\\
&=X_{gh}
\end{align}
$$
Now we can obviously see that since $G$ is a group we can get any vector $X(g)=X_g$
$$
X_g = L_{g,*}X_e
$$
So we can crate a map $ev$ given by
$$
ev:\mathfrak{g}\to T_eG\\
X\in \mathfrak{g}\mapsto X_e
$$
From here on, it is not difficult to show that $ev$ is a vector space isomorphism. (Also one can show that vector space isomorphisms induce diffeomorphisms of smooth manifolds, just saying) $\Box$



## The Exponential Map

Ok now, that we have some very basic constructions defined, it is time to define the exponential map. The expoenntial map is a way to move on the Lie Group by moving on the Lie algebra. Why is this nice? Because it relates the mupltiplication operation of the Lie Group with the addition operation of the Lie algebra as a vector space. This way we can apply a bunch of matrix crap on our Lie algebra and then project them on the lie group. 



**<u>Definition:</u>** Given a vector field $X \in \mathfrak{X}(M)$ we can define an **integral curve** $\gamma: \mathbb{R} \supset I \to M$  through $p\in M$ iff

1. $$\gamma(0) = p$$
2. $$\dot{\gamma}(t) = X_{\gamma(t)} \ \forall t\in I$$

**<u>*Theorem:*</u>** Give a vector field $X$ and a point $p \in M$ we can always find a unique integral curve of $X$ at $p$. 

***<u>Theorem:</u>*** For all $p\in M$ there exists a neigborhood $U \subset M$ of $p$ and an open interval $I$ around $0$ such that the integral curves $\gamma_q$ are defined on $I$ for all $q \in U$. We can even create a map
$$
\begin{align}
\phi_U:U\times I &\to M\\
(q,t) &\mapsto \gamma_q(t)
\end{align}
$$
that is differentiable and is called a **local flow** of $X$.

We can also create a global flow $\phi$ of $X$ on $M$ iff the manifold is closed (compact without boundary). In this case $\phi(\cdot,t):M\to M$ is a diffeomorphism.

***<u>Theorem (Integral Curves of Left-Invariant Vector Fields):</u>*** Consider a lie group $G$ and its Lie algebra $\mathfrak{g}$. Let
$$
\begin{align}
\phi_X:\mathbb{R} \supset I &\to G\\
t &\mapsto \phi_X(t) \in G
\end{align}
$$
to be the maximal integral curve of a vector field $X \in \mathfrak{g}$ that passes through the neutral element $e \in G$. Then the following hold:

1. $\phi_X$ is define on all $\mathbb{R}$.

2. $\phi_X: \mathbb{R} \to G$ is a Lie Group Homomorphism
   $$
   \phi_X(s+t)=\phi_X(s)\cdot\phi_X(t)\ \forall s,t \in \mathbb{R}
   $$

3. $\phi_{sX}(t) = \phi_X(st),\ \forall s,t \in \mathbb{R}$

**<u>Definition:</u>** We define the **exponential map**
$$
\begin{align}
\exp: \mathfrak{g} &\to G\\
X&\mapsto \exp(X) = \phi_X(1)
\end{align}
$$
where $\phi_X$ is the integral curve of some left invariant vector field $X \in \mathfrak{g}$ through the identity $e \in G$.

**<u>*Proposition:*</u>** The exponential map has the following properties

1. $\exp(0) = e$
2. $\exp((s+t)X) = \exp(sX)\cdot \exp(tX)$
3. $\exp(-X) = (\exp X)^{-1}$

***Proof:*** The proof of the statement above is simply through direct calculation. 

1. $\exp(0) = \phi_0(1) = \phi_{0\cdot X}(1) = \phi_{X}(0) = e$ where $X \in \mathfrak{g}$ is some element of the Lie agebra.

2. Follows directly from the fact that the maximal integral curve is a lie group homomorphism. 

3. By using (2) and (1) we have:

   $e=\exp{0}=\exp((1-1)X) = \exp(X)\cdot \exp(-X)$ 

   which implies the statement.

