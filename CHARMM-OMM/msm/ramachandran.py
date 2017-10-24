        phi = md.compute_phi(w)
        psi = md.compute_psi(w)
        ax = ramachandran(phi[1], psi[1], w)
        plt.show()
