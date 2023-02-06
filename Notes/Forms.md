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







