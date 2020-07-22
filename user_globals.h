#ifndef g_kw
    #define g_kw 2000.*CONST_PI  // For changing k
#endif

#ifndef g_m
    #define g_m 0.0   // For changing metallicity slope
#endif

#ifndef g_cz
    #define g_cz 0.6   // For changing metallicity intercept
#endif

#ifndef g_mx
    #define g_mx 0.0   // For changing X slope
#endif

#ifndef g_cx
    #define g_cx 1.01   // For changing X intercept
#endif

#ifndef g_dp
    #define g_dp 0.002   // For changing background pressure, to change temperature
#endif

#ifndef g_drho
    #define g_drho 0.000001    // For changing density fluctuation
#endif

#ifndef g_phi
    #define g_phi (0.0/2.0)*CONST_PI    // For adding phase difference to Z in eigenmode
#endif

#ifndef g_ceil
    #define g_ceil 9.0e8    // To put temperature ceiling
#endif

#ifndef g_noch
    #define g_noch 0.0    // Region of no heating and cooling near boundary
#endif

#ifndef g_iso
    #define g_iso -1    // To choose limit=> 1: Isobaric;   -1: Isochoric
#endif
