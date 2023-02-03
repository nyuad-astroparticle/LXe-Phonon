# Representation Theory

Studying representation theory of groups for particle physics. Math notes

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
3. $[a,[b,c]] + [c,[a,b]] + [b,[c,a]] = 0,\ \forall a,b,c\in \mathbb{V}$