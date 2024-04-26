# Abstract
Homomorphic encryption is a special public-key encryption scheme (PKC), which allows any third party to perform computations on the data while it remains encrypted, despite not having the secret decryption key. Therefore, homomorphic encryption scheme has a wide range of application scenarios, such as protecting the data privacy of users in cloud computing. However, although the efficiency of the homomorphic encryption scheme continues to improve, there are still huge problems, such as ciphertext expansion and time complexity, resulting in its poor efficiency. For example, Cheon et al. proposed CKKS scheme for real approximation based on the problem of learning with error over polynomial rings (RLWE), which can be used in real-world scenarios such as machine learning with privacy protection. However, its plaintext polynomial length is twice the message length. Kim et al. proposed a variation of CKKS scheme based on the conjugate-invariant ring (CIR), which reduces the ratio of plaintext to message. However, this scheme does not support the rotation operation of messages. In view of this, we have realized and improved the scheme proposed by Kim et al. First of all, we adopt the message encoding method proposed by Cheon et al., which makes the scheme support the rotation operation of message slots. Secondly, we propose a transformation method of polynomials over conjugate invariant rings and use fast Fourier transform (FFT) and number theoretic transform (NTT) to accelerate the multiplication of polynomials whose coefficients are real and integer, respectively. Finally, we use C++ language to complete the implementation of the algorithm library of the scheme, including the operations of the underlying polynomials and the cryptography scheme of the upper layer. It supports homomorphic addition, multiplication and rotation operations of the encrypted messages so that it can be applied to more real scenes.


## Backgrounds

### Cyclotomic Ring

To introduce the cyclotomic field and the cyclotomic ring, let us first define the $n$-th cyclotomic polynomial.

**Definition** For any positive integer $n$, the $n$-th cyclotomic polynomial $\phi_n(x)$ is the unique monic polynomial with integer coefficients that divides $x^n - 1$, and for any $k < n$, $x^k - 1$ is not divisible by it.

There is another intuitive definition for cyclotomic polynomials:

**Definition** For any positive integer $n$, the $n$-th cyclotomic polynomial $\phi_n(x)$ is the monic polynomial with integer coefficients, which is the minimal polynomial over the field of numbers for all 𝑛-th primitive roots of unity.

An important property of cyclotomic polynomials is:

**Property** $\prod_{d \mid n} \phi_x(x) = x^n - 1$.

In this post, we use the $M$-th cyclotomic polynomial, where $M = 2^m$. Let $\xi$ be one of the primitive roots of unity; for all elements whose order is 𝑀, they take the form $\xi^1, \xi^3, \ldots, \xi^{M-1}$. According to the second definition of cyclotomic polynomials, we have

$$
\begin{equation}
\phi_M(X) = \prod_{k=1}^{M/2}(X - \xi^{2k - 1}) = X^N + 1
\end{equation}
$$

Here, the degree of $\phi_M(X)$ is $N = \varphi(M) = \frac{M}{2}$, where $\varphi(\cdot)$ is the Euler's totient function.
Since $\phi_M(X)$ is an irreducible polynomial over the real numbers, $\mathcal{F} = \mathbb{Q}[X] / (\phi_M(X))$ forms an extension field of the rational number field $\mathbb{Q}$. We define \mathcal{F} as the cyclotomic field. Meanwhile, $\mathcal{R} = \mathbb{Z}[X] / \phi_M(X)$ is defined as the ring of integers over the field $\mathcal{F}$, which we define as the cyclotomic ring. In this thesis, the entire cryptographic scheme is defined on the residue ring of the cyclotomic ring $\mathcal{R}$ modulo $q$, denoted as $\mathcal{R}_q = \mathcal{R}/q\mathcal{R}$.

### Canonical Embedding

First, we present the definition of embedding in field theory:  

**Definition** An embedding from one field $\mathcal{E}$ to another field $\mathcal{F}$ is a ring homomorphism $\tau: \mathcal{E} \rightarrow \mathcal{F}$. 

Because the kernel of a ring homomorphism is an ideal, and since a field, as a special kind of ring, only has trivial ideals (namely, $\{0\}$ and the field itself), and because $\tau(1) = 1$, it is impossible for this ideal to be the field itself. Consequently, all ring homomorphisms on fields are injective. That is, an embedding is an isomorphism between a field and a subfield of another field. In this thesis, since $M = 2^m$ and $j = 1, 3, \ldots, 2N - 1$, $\xi_j = \xi^j$ are all the roots of $\phi_M(X)$. At this point, for the polynomial $a(X)$ over $\mathcal{F} = \mathbb{Q}[X] / (\phi_M(X))$,
$$\tau_j(X) = \sum_{i = 0}^{N - 1}a_j\xi_j^i$$

Because this mapping is a homomorphic mapping from $\mathcal{F}$ to $\mathbb{C}$, they are all embeddings. The set $\{\tau_j\}_{j = 1, 3, \ldots, 2N−1}$ constitutes a canonical embedding from $\mathcal{F}$ to $\mathbb{C}^{\mathbb{N}}$. We can naturally extend this operation to the field $\mathcal{E} = \mathbb{R}[X] / (X^N + 1)$, obtaining a canonical embedding from $\mathcal{E}$ to $\mathbb{C}^{\mathbb{N}}$. The $l_\infty$ norm of $\tau(a)$ is called the canonical embedding norm of $a$, denoted as $\|a\|^{can}_\infty = \|\tau(a)\|_\infty$. The canonical embedding norm satisfies the following properties:

**Property** For $a, b \in \mathcal{E}$, it holds that $\|a \cdot b\|_{\infty}^{can} \leq \|a\|_{\infty}^{can} \cdot \|b\|_{\infty}^{can}$

**Property** For $a \in \mathcal{E}$, it holds that $\|a\|_{\infty}^{can} \leq \|a\|_1$.

**Property** For a constant $c_M$ related to the ring and valid for all $a \in \mathcal{E}$, it holds that $\|a\|_\infty \leq c_M \cdot \|a\|_{\infty}^{can}$.

The ring constant $c_M = \|V^{-T}\|_\infty$, where $V$ is the Vandermonde matrix formed by the $M$-th principal roots of unity, and the norm of matrix $U$ is defined as $\max_{0 \leq i < N}\{\sum_{i = 0}^{N - 1}|u_{ij}|\}$.

### The Ring of Integers
In field theory, numbers on $\mathbb{Z}$ are referred to as rational integers. Then, the concept of algebraic integers is defined as follows:

**Definition** In the field $\mathcal{F}$, $\alpha$ is an algebraic number if its minimal polynomial has coefficients that are rational integers. If this condition is met, $\alpha$ is called an algebraic integer.

All algebraic integers form a ring, which we call the ring of integers, denoted by $\mathcal{O}_{\mathcal{F}}$. For a set of basis $\{b_0, b_1, \ldots, b_{n - 1}\}$ in $\mathcal{O}_{\mathcal{F}}$, where $n$ is the rank of the number field $\mathcal{F}$, elements in $\mathcal{O}_{\mathcal{F}}$ can be expressed linearly in terms of this basis.

---

## Breif Introduction to CKKS

In this section, we provide a complete introduction to the RNS version of the CKKS scheme.
- $CKKS.KeyGen(1^\lambda)$: First, select $L$ as the number of prime numbers and the depth of homomorphic multiplication, and choose $L$ primes $q_0, q_1, \ldots, q_{L-1}$, satisfying $q_i \in \left(2^p - \frac{1}{2^\eta}, 2^p + \frac{1}{2^\eta}\right)$. Next, based on the security parameter $\lambda$, choose the polynomial degree $N = 2^n$. Finally, generate the secret $s \leftarrow \mathcal{HWT}$, set the secret key as $\text{sk} = (1, s)$, and compute the public key and evaluation key as $\text{pk} = (-a_{\text{pk}} s + e_{\text{pk}}, a_{\text{pk}})$, $\text{evk} = (-a_{\text{evk}} s + e_{\text{evk}} + P s^2, a_{\text{evk}})$, where $a_{\text{pk}} \leftarrow \mathcal{U}(\mathcal{R}_q)$, $a_{\text{evk}} \leftarrow \mathcal{U}(\mathcal{R}_q)$, $e_{\text{pk}}, e_{\text{evk}} \leftarrow \mathcal{DG}(\sigma)$, $P = p_0 p_1 \ldots p_{K-1}$ is a square-free integer. Output $(\text{pk}, \text{evk}, \text{sk})$.
- $CKKS.Encode(\mathbf{v}, p)$: Input a message of length at most $\frac{N}{2}$ and a scale factor $\text{scale} = 2^p$, first extend it to an element in the linear space $\mathcal{H}$, then use the inverse operation of canonical embedding to obtain a real polynomial and scale it by $\text{scale}$, finally round to get an integer-coefficient polynomial and output it, i.e., $m = \left\lfloor 2^p \tau^{-1} (\pi^{-1} (\mathbf{v})) \right\rceil_{\mathcal{R}}$.
- $CKKS.Decode(m, p)$: Input a plaintext and a scale, use canonical embedding to first convert it into a message vector and scale it, finally output the first half, i.e., $\mathbf{v} = \pi \left( \frac{1}{2^p} \tau (m) \right)$.
- $CKKS.Enc_{\mathtt{pk}}(m)$: Input a plaintext, first generate a random number $r \leftarrow ZO$, noise $e_0, e_1 \leftarrow \mathcal{DG}(\sigma)$, then output the ciphertext $\text{ct} = r \cdot \text{pk} + (m + e_0, e_1)$.
- $CKKS.Dec_{\mathtt{sk}}(\text{ct})$: Input a ciphertext, output the decrypted plaintext, i.e., $m = \langle \text{ct}, \text{sk} \rangle$.
- $CKKS.Add(\text{ct}_1, \text{ct}_2)$: Input two ciphertexts, output the addition ciphertext as $\text{ct}_{\text{add}} = \text{ct}_1 + \text{ct}_2$.
- $CKKS.Mult_{\mathtt{evk}}(\text{ct}_1, \text{ct}_2)$: Input two ciphertexts, denoted $\text{ct}_1 = (b_1, a_1)$, $\text{ct}_2 = (b_2, a_2)$, let $d_0 = b_1b_2$, $d_1 = a_1b_2 + a_2b_1$ and $d_2 = a_1a_2$. First, use the Chinese Remainder Theorem to merge and decompose $d_2$ into $d_2' \in \mathcal{R}_{P_q}$, then output the multiplication ciphertext as $\text{ct}_{\text{mult}} = (d_0, d_1) + \left\lfloor P^{-1} d_2' \cdot \text{evk} \right\rceil_{\mathcal{R}_q}$.
- $CKKS.Rescale(\text{ct}, l')$: Input a ciphertext needing rescaling and the rescaling layer, where the original ciphertext is an element over $\mathcal{R}_{q_0 \ldots q_l}$, output the rescaled ciphertext $\text{ct}' \leftarrow \left\lfloor 2^{-l'} \cdot \text{ct} \right\rceil_{\mathcal{R}_{q_0 \ldots q_{l-l'}}}$.

---

## Our Algorithm

In this chapter, we introduce improvements to the CKKS scheme. Unlike the original CKKS scheme, the message vector space for this scheme is $\mathbb{R}^{N/2}$. First, in most real-world scenarios, only the handling of real numbers is required. Secondly, even when dealing with complex numbers, we only need to treat the real and imaginary parts as two independent real numbers and encrypt them separately, then calculate according to the rules of complex number operations. When the message space is restricted to real number vectors, because the conjugate of a real number is itself, the part of the message expansion is the same as the original vector. Intuitively, we only need polynomials of length $N/2$, thus achieving higher efficiency. However, as Theorem 3.1 proved, only when the message is in the space $\mathcal{H}$ can we ensure the polynomial coefficients are real. Therefore, we need to use a new polynomial ring structure. Kim et al. [2] proposed building CKKS on the Conjugation-invariant ring (CIR), changing the polynomial structure from $\mathbb{Z}_q[X] / (X^N + 1)$ to $\mathbb{Z}_q[X] / (X^{N/2} + X^{-N/2})$, which allows the CKKS to shorten the polynomial length when only dealing with real numbers. In this paper, we further extend Kim et al.'s scheme in two aspects. On one hand, we further change the form of canonical embedding, thus implementing homomorphic rotation operations. This operation can move the message vector any number of positions to the left or right, equivalent to completing a cyclic shift. This operation plays an important role in implementing matrix multiplication and many other applications. On the other hand, we designed a form transformation of polynomials on the CIR, enabling polynomial calculations in this ring space to be accelerated by the NTT algorithm.

### New Polynomial Representation and Canonical Embedding Operations

In this section, we introduce the representation of polynomial elements in the CIR and the corresponding canonical embedding operations.

In the original CKKS scheme, when encoding real number message vectors, the first step requires expanding them from $\mathbb{R}^{N/2}$ to $\mathbb{R}^N$ such that for $i = 0, \ldots, N/2 - 1, v_{i+N/2} = \overline{v_i} = v_i$. This requires that the polynomial $a(X)$ when substituted into $\xi_{2k+1}$ and $\xi_{M-2k-1} = \xi^{-1}_{2k+1}$, yields the same results. That is, the plaintext space becomes a subspace of $\mathcal{R}$ $\{ a \in \mathcal{R} | a(X) = a(X^{-1}) \}$. Due to the characteristics of polynomial operations, it is easy to see this is a subring of $\mathcal{R}$, which we call the CIR. To facilitate polynomial calculations, we need a simple form to represent all elements on the CIR.

We define the number field $\mathcal{F}$ as $\{ a \in \mathbb{Q}[X] / (X^N + 1) | a(X) = a(X^{-1}) \}$. Additionally, based on the knowledge of algebraic structures, we know when both $(X^{N/2} + X^{-N/2})$ and $(X^N + 1)$ are irreducible polynomials over $\mathbb{Q}$ and $\zeta, \xi$ are respectively roots of $(X^{N/2} + X^{-N/2})$ and $(X^N + 1)$, we have $\mathbb{Q}(\zeta) = \mathbb{Q}[X] / (X^{N/2} + X^{-N/2})$ and $\mathbb{Q}(\xi) = \mathbb{Q}[X] / (X^N + 1)$. Below we present a theorem:

**Theorem** The number field $\mathbb{Q}[X] / (X^{N/2} + X^{-N/2}) = \mathcal{F}$.

**Proof** Because $\mathbb{Q}(\zeta) = \mathbb{Q}[X] / (X^{N/2} + X^{-N/2})$, we only need to prove $\mathbb{Q}(\zeta) = \mathcal{F}$. It is easy to see that $\mathbb{Q}(\zeta) \subseteq \mathcal{F} \subseteq \mathbb{Q}(\xi)$; next, we only need to prove $\mathcal{F} \subseteq \mathbb{Q}(\zeta)$. We define the cyclic group $G = \{\text{id}, \tau^{-1}\}$, where the generator $\tau^{-1} : \xi \to \xi^{-1}$, then $\mathcal{F}$ is the fixed point of $G$ acting on $\mathbb{Q}(\xi)$, and $


Because $\mathbb{Q}(\zeta) = \mathbb{Q}[X] / (X^{N/2} + X^{-N/2})$, we only need to prove $\mathbb{Q}(\zeta) = \mathcal{F}$. It is easy to find that $\mathbb{Q}(\zeta) \subseteq \mathcal{F} \subseteq \mathbb{Q}(\xi)$; next, we only need to prove $\mathcal{F} \subseteq \mathbb{Q}(\zeta)$. We define the cyclic group $G = \{\text{id}, \tau_{-1}\}$, where the generator $\tau_{-1} : \xi \to \xi^{-1}$, then $\mathcal{F}$ is the fixed point of $G$ acting on $\mathbb{Q}(\xi)$, and $[\mathbb{Q}(\xi) : \mathcal{F}] = |G| = 2$. Thus, it is equivalent to proving $[\mathbb{Q}(\xi) : \mathbb{Q}(\zeta)] \leq 2$. Clearly, both $\xi$ and $\xi^{-1}$ are roots of $(X^{N/2} + X^{-N/2})$, so $\xi + \xi^{-1} \in \mathbb{Q}(\zeta)$. For a polynomial $X^2 - (\xi + \xi^{-1}) \cdot X + 1$ with coefficients in $\mathbb{Q}(\zeta)$, it is found that $\xi$ is a root of it. Therefore, the minimal polynomial of $\xi$ in the polynomial ring $\mathbb{Q}(\zeta)[X]$ has degree at most 2, i.e., $[\mathbb{Q}(\xi) : \mathbb{Q}(\zeta)] \leq 2$. Therefore, $\mathcal{F} \subseteq \mathbb{Q}(\zeta)$. Combining these points, $\mathbb{Q}(\zeta) = \mathcal{F}$.

Since the polynomials over $\mathbb{Q}(\zeta) = \mathbb{Q}[X] / (X^{N/2} + X^{-N/2})$ can be expressed in the following form:

$$
a(X) = a_0 + \sum_{i=1}^{N/2-1} a_i (X^i + X^{-i})
$$

Due to the equivalence of $\mathbb{Q}(\zeta)$ and $\mathcal{F}$, we obtain the representation on $\mathcal{F}$. It is important to note that the CIR ring is the ring of integers over $\mathcal{F}$, so elements of the CIR ring can also be represented in this form. In subsequent chapters, we use $\mathcal{R}'$ to denote the CIR ring, $\mathcal{F}'$ to denote $\mathbb{Q}[X] / (X^{N/2} + X^{-N/2})$, and $\mathcal{E}'$ to denote $\mathbb{R}[X] / (X^{N/2} + X^{-N/2})$.

### New Canonical Embedding Operations and Rotation Operations

To encode a message vector onto a polynomial in the CIR ring, it is necessary to substitute the roots of the irreducible polynomial $X^{N/2} + X^{-N/2}$, which forms a basis for canonical embedding. Given that
$$x^{N/2} + x^{-N/2} = x^{-N/2}(x^N + 1)$$

thus $\zeta_1, \zeta_3, \dots, \zeta_{2N-1}$ are roots of this equation. However, for any $1 \leq i \leq 2N - 1$ and polynomial $a \in \mathcal{R}$, $a(\zeta_i) = a(\zeta_{2N-i})$. Moreover, for $N = 2^n > 0$ and all odd $i$, when $i = 1 \mod 4$, $2N - i = 3 \mod 4$, and when $i = 3 \mod 4$, $2N - i = 1 \mod 4$. Thus, $(a(\zeta_{4j+1}))_{j=0}^{N/2-1}$ is an injection from a polynomial to $\mathbb{R}^{N/2}$. The new canonical mapping is defined as:
$$\tau(a) = (\tau_1(a), \tau_5(a), \dots, \tau_{2N-3}(a))$$

For the homomorphic rotation operations on message vectors, Cheon et al. mentioned a rotation method. However, the aforementioned canonical embedding is not compatible with it. To apply it to $\mathcal{R}'$, combined with the techniques from Cheon et al.'s paper, we propose a variant of the canonical embedding on $\mathcal{R}'$. A fact: $5j \mod M$ can traverse all $4k + 1 \mod M$ terms, hence we use $5j$ instead of $(4k + 1)$, defining the variant of the canonical embedding operation as:
$$\tau(a) = (a(\zeta^{5^0}), a(\zeta^{5^1}), a(\zeta^{5^{N/2-1}}))$$

Additionally, assume we can construct a mapping $\kappa_k : a(X) \rightarrow a(X^k)$. Using the mapping $\kappa_{5r}$, we can construct an operation that rotates a message vector by $r$ positions, as follows:
$$a'(\zeta^{5j}) = a(\zeta{5^j \cdot 5^r}) = a(\zeta^{5^{j+r}})$$

Below, we present the construction of $\kappa_k(\cdot)$:
$$a(X^k) = a_0 + \sum_{i=1}^{N/2} a_i (X^{ik} + X^{-ik})$$

Let $\kappa_k(a) = a'$, where for $a'(X) = a'_0 + \sum_{i=1}^{N/2-1} a'_i (X^i + X^{-i})$, the coefficients of $a'(X)$ relate to the original polynomial coefficients as $a'_{ik} = a_i$. Note, the correctness of this construction requires $k$ and $M$ to be coprime.

For plaintext, the rotation operation is simply performing the aforementioned rotation transformation directly on the polynomial. For ciphertext, it is given by:
$$
\begin{align}
\kappa(m) 
&\approx \langle \kappa(\text{ct}), \kappa(\text{sk}) \rangle \\
&= \kappa(b) + \kappa(a) \kappa(s)
\end{align}
$$

Thus, similar to ciphertext multiplication, we need to perform a key-switching operation. Let the rotation key be $\text{rtk} = (b_{\text{rtk}} = -a_{\text{rtk}}s + e_{\text{rtk}} + P \cdot \kappa(s), a_{\text{rtk}})$, where $a_{\text{rtk}} \leftarrow \mathcal{U}(\mathcal{R}'_{P_q})$, $e_{\text{rtk}} \leftarrow \mathcal{DG}(\sigma)$, then we define the rotated ciphertext $\text{ct}_{\text{rot}} = (\kappa(b), 0) + P^{-1} \kappa(a) \cdot \text{rtk}$. The ciphertext decrypts to:

$$
\begin{align}
    \langle \text{ct}_{\text{rot}}, \text{sk} \rangle \\
    &= \kappa(b) + P^{-1}\kappa(a) \cdot (b_{\text{rtk}} + a_{\text{rtk}}s) \\
    &= \kappa(b) + P^{-1}\kappa(a) \cdot (P \kappa(s) + e_{\text{rtk}}) \\
    &= \kappa(b) + \kappa(a) \kappa(s) + P^{-1}\kappa(a) e_{\text{rtk}} \\
    &= \kappa(m) + \kappa(e) + P^{-1}\kappa(a) e_{\text{rtk}}
\end{align}
$$
From this, we obtain the plaintext after rotation, where the noise generated by the rotation operation is $e_{\text{rot}} = \kappa(e) + P^{-1}\kappa(a) e_{\text{rtk}}$.

### Fast Fourier Transform and Fast Number Theoretic Transform

In this section, we introduce the form transformation of polynomials for the NTT algorithm. The bottleneck in polynomial operations is polynomial multiplication. If the degrees of polynomials $a$ and $b$ are $n$ and $m$ respectively, multiplying them using the naive method would have a time complexity of $O(nm)$. Such complexity is intolerable when the polynomial degrees are very high. A commonly used efficient algorithm is the Fast Fourier Transform (FFT). With the FFT algorithm, the complexity of polynomial multiplication can be reduced to $O(n \log n)$. The general principle is as follows:

Firstly, we define two types of polynomial representations: The first is the coefficient representation, where for a polynomial $a(X) = \sum_{i=0}^{n} a_iX^i$, we represent $a(X)$ by $a_c = (a_0, a_1, \ldots, a_n)$. The other is called the point-value representation, where for a polynomial $a(X)$, we take a set of primitive roots of unity $\zeta_0, \zeta_1, \ldots, \zeta_n \in \mathbb{C}$, and represent $a(X)$ as $a_v = (a(\zeta_0), a(\zeta_1), \ldots, a(\zeta_n))$. The following property holds for the point-value representation of polynomials:

**Property** For polynomials $a(X)$ and $b(X)$, we have:
$$(ab)_v = (a(\zeta_0)b(\zeta_0), a(\zeta_1)b(\zeta_1), \ldots, a(\zeta_n)b(\zeta_n))$$

The FFT algorithm is a fast method for converting a polynomial from coefficient representation to point-value representation, with a complexity of $O(n \log n)$, while the inverse FFT algorithm is the reverse process of the FFT. Correspondingly, for polynomials over $\mathbb{Z}_q[X] / (X^N + 1)$, when $q$ is a prime, its primitive roots also have properties similar to the primitive roots of unity in the complex field. Therefore, the variant of FFT in finite polynomial rings is the Number Theoretic Transform (NTT).

However, due to the unique polynomial representation in the CIR ring $\mathcal{R}'$, the FFT and NTT algorithms cannot be directly applied.

Thus, we need to transform its form as follows:

$$
\begin{align}
a(\zeta^{5^j}) 
&= a_0 + \sum_{i=1}^{N/2-1} a_i (\zeta^{i \cdot 5^j} + \zeta^{-i \cdot 5^j}) \\
&= a_0 + \sum_{i=1}^{N/2-1} a_i \zeta^{i \cdot 5^j} + \sum_{i=1}^{N/2-1} a_i \zeta^{-i \cdot 5^j} \\
&= a_0 + \sum_{i=1}^{N/2-1} a_i \zeta^{i \cdot 5^j} + \sum_{i=1}^{N/2-1} a_{N/2-i} \zeta^{(i-N/2) \cdot 5^j} \\
&= a_0 + \sum_{i=1}^{N/2-1} \left( a_i \zeta^{i \cdot 5^j} + a_{N/2-i} \zeta^{(i-N/2) \cdot 5^j} \right) \\ 
&= a_0 + \sum_{i=1}^{N/2-1} \left( a_i + a_{N/2-i} \zeta^{-N/2 \cdot 5^j} \right) \zeta^{i \cdot 5^j} \\
&= a_0 + \sum_{i=1}^{N/2-1} (a_i + a_{N/2-i} \zeta^{-N/2}) \zeta^{i \cdot 5^j} \\
&= a'(\zeta^{5^j})
\end{align}
$$

Where $a'(X) = (a_0, a_1 + a_{N/2-1} \zeta^{-N/2}, \ldots, a_{N/2-1} + a_1 \zeta^{-N/2})$ can be used in the FFT and NTT algorithms. Thus, when computing the point values of $a(X)$, we first need to transform it into $a'(X)$, then apply the FFT (or NTT) algorithm.

When converting a polynomial from point-value representation to coefficient representation, we first call the FFT algorithm (or NTT algorithm) to obtain the polynomial $a'(X)$, and then based on the formula $a(X) = (a'_0, a'_1 + a'_{N/2-1} \zeta^{N/2}, \ldots, a'_{N/2-1} + a'_1 \zeta^{N/2})$ to get the polynomial $a(X)$ on the CIR ring.

## Implementation
We have implemented the improvements to the RNS-CKKS scheme using the C++ programming language, naming it [SOPHIE](https://github.com/dreamATD/SOPHIE) (ShOrt Polynomial HomomorphIc Encryption). During the development of this cryptographic library, we used GitHub for code version management and the Catch library for unit testing. To enhance the efficiency of our implementation, we utilized the modern language feature of move semantics introduced in C++11.

# References
1. [Homomorphic Encryption for Arithmetic of Approximate Numbers](https://eprint.iacr.org/2016/421)

2. [Approximate Homomorphic Encryption over the Conjugate-Invariant Ring](https://eprint.iacr.org/2018/952)

3. [Bootstrapping for Approximate Homomorphic Encryption](https://eprint.iacr.org/2018/15)
