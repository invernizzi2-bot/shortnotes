#!/usr/bin/env python3

import math


# ============================================================
# Ellissoide internazionale 1924 (Hayford)
# ============================================================

a = 6378388.0          # semiasse maggiore [m]
inv_f = 297.0          # inverso dello schiacciamento

f = 1.0 / inv_f
e2 = 2.0 * f - f * f
e = math.sqrt(e2)


# ============================================================
# Coefficienti della proiezione Equal Earth
# ============================================================

A1 = 1.340264
A2 = -0.081106
A3 = 0.000893
A4 = 0.003796

M = math.sqrt(3.0) / 2.0


# ============================================================
# Conversione da latitudine geodetica a latitudine authalica
# ============================================================

def geodetic_to_authalic(phi):
    """
    Converte la latitudine geodetica phi [radianti]
    nella latitudine authalica beta [radianti].
    """

    sin_phi = math.sin(phi)

    q = (1.0 - e2) * (
        sin_phi / (1.0 - e2 * sin_phi * sin_phi)
        - (1.0 / (2.0 * e))
        * math.log(
            (1.0 - e * sin_phi) /
            (1.0 + e * sin_phi)
        )
    )

    # Valore di q al polo
    qp = (
        1.0
        - (1.0 - e2) / (2.0 * e)
        * math.log(
            (1.0 - e) /
            (1.0 + e)
        )
    )

    # Relazione tra q e latitudine authalica
    sin_beta = q / qp

    # Protezione contro piccoli errori numerici
    sin_beta = max(-1.0, min(1.0, sin_beta))

    return math.asin(sin_beta)


# ============================================================
# Proiezione Equal Earth
# ============================================================

def equal_earth(lam_deg, phi_deg, lon0_deg=0.0):
    """
    Calcola le coordinate x, y della proiezione Equal Earth.

    Parametri
    ---------
    lam_deg : float
        Longitudine geodetica [gradi].

    phi_deg : float
        Latitudine geodetica [gradi].

    lon0_deg : float
        Longitudine del meridiano centrale [gradi].

    Restituisce
    -----------
    x : float
        Coordinata Est [m].

    y : float
        Coordinata Nord [m].
    """

    # Conversione da gradi a radianti
    lam = math.radians(lam_deg - lon0_deg)
    phi = math.radians(phi_deg)

    # Latitudine geodetica -> latitudine authalica
    beta = geodetic_to_authalic(phi)

    # Calcolo dell'angolo theta
    sin_theta = M * math.sin(beta)

    # Protezione contro errori numerici
    sin_theta = max(-1.0, min(1.0, sin_theta))

    theta = math.asin(sin_theta)

    # Polinomio della Equal Earth
    F_theta = (
        A1 * theta
        + A2 * theta**3
        + A3 * theta**7
        + A4 * theta**9
    )

    # Derivata del polinomio
    D_theta = (
        A1
        + 3.0 * A2 * theta**2
        + 7.0 * A3 * theta**6
        + 9.0 * A4 * theta**8
    )

    # Coordinate Equal Earth
    x = (
        a
        * lam
        * math.cos(theta)
        / (M * D_theta)
    )

    y = a * F_theta

    return x, y


# ============================================================
# Esempio di utilizzo
# ============================================================

if __name__ == "__main__":

    lam = float(input("Longitudine λ [gradi]: "))
    phi = float(input("Latitudine φ [gradi]: "))

    if not -180.0 <= lam <= 180.0:
        raise ValueError(
            "La longitudine deve essere compresa "
            "tra -180° e +180°."
        )

    if not -90.0 <= phi <= 90.0:
        raise ValueError(
            "La latitudine deve essere compresa "
            "tra -90° e +90°."
        )

    x, y = equal_earth(lam, phi)

    print()
    print(f"x = {x:.6f} m")
    print(f"y = {y:.6f} m")
