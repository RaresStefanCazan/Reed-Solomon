import sympy
import random
from itertools import combinations
import sys
import time
sys.stdout.reconfigure(encoding='utf-8')


def generate_large_numbers(num_bits=161):
    lower_bound = 1 << (num_bits - 1)
    upper_bound = (1 << num_bits) - 1
    return sympy.randprime(lower_bound, upper_bound)

def text_to_ascii_list(message):
    #transform mesajul in lista de coduri ASCII
    return [ord(c) for c in message]

def split_list(lst, block_size=20):
    #Impart lista de coduri ASCII in blocuri de 20
    return [lst[i:i+block_size] for i in range(0, len(lst), block_size)]

def block_list_to_int(block):
    # 256^20 = 2^160 deci toate vor fi mai mici ca p
    result = 0
    for code in block:
        result = result * 256 + code
    return result

def polynom(coeffs, x, p):
    #horner
    val = 0
    power = 1
    for a in coeffs:
        power = (power * x) % p
        val = (val + a * power) % p
    return val

def encode(message, p, block_size=20):
    #codificarea

    ascii_list = text_to_ascii_list(message)
    blocks = split_list(ascii_list, block_size)
    coeffs = [block_list_to_int(block) for block in blocks]
    k = len(coeffs)
    n = k + 2
    codeword = [polynom(coeffs, x, p) for x in range(1, n+1)]
    return codeword, k

#incep calcul coef liber

def modinv(a, p):
    #Calculeaza inversul modular cu fermat a^p-2 = a^-1 % p
    return pow(a, p-2, p)

def free_coef_naive(z, A, p):
    #formula din pdf

    fc = 0
    for i in A:
        zi = z[i-1]
        product = 1
        for j in A:
            if j != i:
                product = (product * j) % p
                diff = (i - j) % p
                inv_diff = modinv(diff, p)
                product = (product * inv_diff) % p
        fc = (fc + zi * product) % p
    return fc

def free_coef_k(z, A, p):
    fc = 0
    for i in A:
        zi = z[i-1]
        numa_prod = 1
        numi_prod = 1
        for j in A:
            if j != i:
                numa_prod = (numa_prod * j) % p
                numi_prod = (numi_prod * ((j - i) % p)) % p
        # fac o singura inversare la final pentru tot produsul de la numitor
        inv_den = modinv(numi_prod, p)
        prod_i = (numa_prod * inv_den) % p
        fc = (fc + zi * prod_i) % p
    return fc


def free_coef_one_inversion(z, A, p):
    #o singura inversare ca mai intai rezolv matematic separat numaratorul si numitorul

    total_num, total_den = 0, 1
    for i in A:
        num_i, den_i = fraction_for_index(i, z, A, p)
        total_num, total_den = add_fractions(total_num, total_den, num_i, den_i, p)
    inv_total_den = modinv(total_den, p)
    return (total_num * inv_total_den) % p

def fraction_for_index(i, z, A, p):
    num = z[i-1]
    den = 1
    for j in A:
        if j != i:
            num = (num * j) % p
            den = (den * ((j - i) % p)) % p
    return num, den

def add_fractions(n1, d1, n2, d2, p):
    n = (n1 * d2 + n2 * d1) % p
    d = (d1 * d2) % p
    return n, d

def lagrange_eval(points, x, p):
   #din pdf evalueaza la pct x

    total = 0
    for j, (xj, yj) in enumerate(points):
        
        num = 1
        den = 1
        for m, (xm, ym) in enumerate(points):
            if m != j:
                num = (num * (x - xm)) % p
                den = (den * (xj - xm)) % p
        #Lj e produsul din formula
        Lj = (num * modinv(den, p)) % p

        
        term = (yj * Lj) % p

        #fac si suma din formula
        total = (total + term) % p

    return total

def correct_single_error_using_freecoef(z, p, A):
  
    n = len(z)
    for e in A:
        # iau toate nodurile mai putin unu
        A_cand = [i for i in A if i != e]
        
        fc_cand = free_coef_naive(z, A_cand, p)
        if fc_cand % p == 0:
            # Daca free_coef este 0, punctele rămase sunt corecte.
            # Calculez valoarea corecta la x = e prin interpolare Lagrange.
            points = [(i, z[i-1]) for i in A if i != e]
            correct_val = lagrange_eval(points, e, p)
            old_val = z[e-1]
            z[e-1] = correct_val
            print(f"Eroare detectată la indexul {e}. Valoare veche = {old_val}, valoare nouă = {correct_val}")
            return z, e
    print("Nicio eroare nu a fost detectată.")
    return z, None


def compare_free_coef_methods(z, A, p, iterations=1000):
    
    t_naive = 0.0
    t_k = 0.0
    t_one_inv = 0.0
    result_naive = None
    result_k = None
    result_one_inv = None

    for _ in range(iterations):
        start = time.perf_counter()
        result_naive = free_coef_naive(z, A, p)
        t_naive += time.perf_counter() - start

        start = time.perf_counter()
        result_k = free_coef_k(z, A, p)
        t_k += time.perf_counter() - start

        start = time.perf_counter()
        result_one_inv = free_coef_one_inversion(z, A, p)
        t_one_inv += time.perf_counter() - start

    print("Rezultate:")    
    print("free_coef_naive =", result_naive)
    print("free_coef_k =", result_k)
    print("free_coef_one_inversion =", result_one_inv)
    print("\nTimp mediu (peste {} execuții):".format(iterations))
    print("free_coef_naive: {:.6f} sec".format(t_naive/iterations))
    print("free_coef_k: {:.6f} sec".format(t_k/iterations))
    print("free_coef_one_inversion: {:.6f} sec".format(t_one_inv/iterations))

def poly_add(A, B, p):
    
    max_len = max(len(A), len(B))
    result = [0]*max_len #lista de 0
    
    for i in range(max_len):
        #daca A are coeficient la i
        ai = A[i] if i < len(A) else 0
        bi = B[i] if i < len(B) else 0
        result[i] = (ai + bi) % p
    
    return result


def poly_mul(A, B, p):
    
    result = [0]*(len(A)+len(B)-1) #lista de 0

    for i in range(len(A)):
        for j in range(len(B)):
            #indexul corespunde gradului polinomului
            result[i+j] = (result[i+j] + A[i]*B[j]) % p

    return result

def poly_scale(A, scalar, p):
    
    return [(a * scalar) % p for a in A]

def reconstruct_polynomial(points, p):

    P = [0]  
    
    for i, (x_i, y_i) in enumerate(points):
        
        Li = [1]  # polinom constant 1
        den = 1
        for j, (x_j, y_j) in enumerate(points):
            if j != i:
                #numarator
                Li = poly_mul(Li, [(-x_j) % p, 1], p)
               #numi
                den = (den * ((x_i - x_j) % p)) % p
        
        #inm cu inversa
        inv_den = modinv(den, p)
        Li = poly_scale(Li, inv_den, p)
        
        #si cu y_i
        Li = poly_scale(Li, y_i, p)
        
        #formez polinomul
        P = poly_add(P, Li, p)
    
    return P

def int_to_block_list(num):
    # transform un intreg intr o lista de bytes
    block = []
    while num > 0:
        block.append(num % 256)
        num //= 256
    return block[::-1]  # inversez ca sa fie ordinea corecta




if __name__ == "__main__":
    p = generate_large_numbers(161)
    print(p)
    message = "Acesta este un mesaj de test pentru a demonstra blocarea in grupuri de 20 de caractere. Pot fi adaugate oricat de multe caractere, dar dupa transformare ascii vor fi impartite in coduri de 20."
    
    # Encoding
    vector_encoded, k = encode(message, p, block_size=20)
    print("Vectorul de codificare:")
    print(vector_encoded)
    
    # Alegem A = {1, 2, ..., n}
    n = len(vector_encoded)
    A = list(range(1, n+1))
    
    # Testez ca fara erori coef liber e 0
    fc_naive = free_coef_naive(vector_encoded, A, p)
    fc_k = free_coef_k(vector_encoded, A, p)
    fc_one_inv = free_coef_one_inversion(vector_encoded, A, p)
    
    print("\nCoeficientul liber (fără eroare):")
    print("free_coef_naive =", fc_naive)
    print("free_coef_k =", fc_k)
    print("free_coef_one_inversion =", fc_one_inv)
    
    #Adaug eroarea
    vector_with_error = vector_encoded.copy()
    error_index = 2  
    bit_to_flip = 3  
    #cu xor
    vector_with_error[error_index] ^= (1 << bit_to_flip)
    
    print("\nVectorul de codificare cu eroare (bit flip la index {}):".format(error_index))
    print(vector_with_error)
    
    # Calculez coeficientul liber pentru vectorul cu eroare: ar trebui sa fie != 0
    fc_err_naive = free_coef_naive(vector_with_error, A, p)
    fc_err_k = free_coef_k(vector_with_error, A, p)
    fc_err_one_inv = free_coef_one_inversion(vector_with_error, A, p)
    
    print("\nCoeficientul liber calculat (cu eroare):")
    print("free_coef_naive =", fc_err_naive)
    print("free_coef_k =", fc_err_k)
    print("free_coef_one_inversion =", fc_err_one_inv)

#corectez eroarea
vector_corrected, error_index = correct_single_error_using_freecoef(vector_with_error, p, A)

print("\nVectorul corectat:")
print(vector_corrected)

if vector_corrected == vector_encoded:
    print("\nCorectare reușită! Vectorul corectat corespunde vectorului original.")
else:
    print("\nCorectare eșuată! Vectorul corectat nu coincide cu cel original.")

fc_after = free_coef_naive(vector_corrected, A, p)
print("\nCoeficientul liber după corectare =", fc_after)


print("\nCompararea timpilor de execuție pentru metodele free_coef:")
compare_free_coef_methods(vector_encoded, A, p, iterations=100)


    #salvez (index,valoare)
points_corrected = [(i, vector_corrected[i-1]) for i in A]
poly_coeffs = reconstruct_polynomial(points_corrected, p)

    # extrag coeficientii a1,..,ak fara a0
message_coeffs = poly_coeffs[1:k+1]

all_bytes = []
for a_i in message_coeffs:
        block_bytes = int_to_block_list(a_i)
        all_bytes.extend(block_bytes)

decoded_message = ''.join(chr(b) for b in all_bytes)
print("Mesaj decodificat:", decoded_message)