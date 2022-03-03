module modamerica

  implicit none

  ! The objective is to draw a map of America in which the countries
  ! appear with areas that are proportional to a target value (real
  ! area, population, PIB, etc). The unknowns are 132 points in R^2,
  ! which define the boundaries of 17 "countries". Each point is
  ! assigned to the boundary of one or more country. The computed area
  ! of each country is calculated as a function of its boundary points
  ! using Green's formula. The constraints of the problem are:
  !
  !               computed area = (scaled) target value
  !
  ! for each country. As initial approximation we took the coordinates
  ! of the 132 points in the New York Times map of America which, of
  ! course, do not fit with target values. The objective function is
  !
  ! 1/2 \sum_{j=1}^{132} \|P_j - Q_j\|_2^2,
  !
  ! where P_1,...,P_132 are the unknowns and Q_1,...Q_132 are the
  ! initial approximation.
  !
  ! Considered "countries" are:
  !
  !  1 Cuba
  !  2 Canada
  !  3 USA (excluding Alaska that is considered in separate)
  !  4 Mexico
  !  5 America Central (Guatemala, Hunduras, Nicaragua, El Salvador, 
  !                     Costa Rica, Panama, and Belize)
  !  6 Colombia
  !  7 Venezuela
  !  8 Guianas (Guyana, Suriname, and French Guiana)
  !  9 Brasil
  ! 10 Ecuador
  ! 11 Peru
  ! 12 Bolivia
  ! 13 Paraguay
  ! 14 Chile
  ! 15 Argentina
  ! 16 Alaska
  ! 17 Uruguay

  ! Number of states and number of points

  integer :: nstates,np
  parameter ( nstates = 17 )
  parameter ( np = 132 )

  ! Number of boundary points for each state

  integer :: nbord(nstates)
  data nbord(1:nstates) / 4, 21, 22, 26, 8, 9, 9, 5, 22, 4, 10, 7, 6, 16, 20, 8, 5 /

  ! Indices of the boundary points per state

  integer pbord(np,nstates)
  data pbord(1: 4, 1) /  49, 50, 51, 52 /, &
       pbord(1:21, 2) /  83, 35, 64, 63, 84, 36, 99,100, 37, 38, 39,101,102, &
                         40, 41, 42,127, 81, 82, 43, 46 /, &
       pbord(1:22, 3) /  33,108,109,110,111,112, 34, 94, 62, 98, 61, 60, 96, &
                         97, 36, 84, 63, 64, 35,120, 56,121 /, &
       pbord(1:26, 4) /  86,113,114, 31, 32, 87,115, 88,103, 89, 85, 34,112, &
                        111,110,109,108, 33,118,104,116,119,105,106,117,107 /, &
       pbord(1: 8, 5) /  30, 29, 91, 90, 32, 31, 92, 93 /, &
       pbord(1: 9, 6) /  28, 25, 23, 22, 55, 95, 24, 29, 30 /, &
       pbord(1: 9, 7) /  22, 21, 20, 19,122,123, 24, 95, 55 /, &
       pbord(1: 5, 8) /  20,126, 18, 17, 19 /, &
       pbord(1:22, 9) /  53, 16, 15, 10,  8,  7,  6, 65,  5, 68, 76, 48,129, &
                         47,125, 17, 18,126, 20, 21, 22, 23 /, &
       pbord(1: 4,10) /  27, 26, 25, 28 /, &
       pbord(1:10,11) / 128, 59, 12, 13, 16, 53, 23, 25, 26, 27 /, &
       pbord(1: 7,12) /  14,  9, 10, 15, 16, 13, 12 /, &
       pbord(1: 6,13) /  54,  7,  8, 10,  9, 70 /, &
       pbord(1:16,14) /   1, 11,132,131, 14, 12, 59,130, 57, 58,  2, 73, 77, &
                         71, 75, 74 /, &
       pbord(1:20,15) /   1,  2, 73, 71, 72, 75, 74, 69, 67, 66,  4,  6,  7, &
                         54, 70,  9, 14,131,132, 11 /, &
       pbord(1: 8,16) /  43, 80, 44, 78,124, 45, 79, 46 /, &
       pbord(1: 5,17) /   4,  3,  5, 65,  6 /

  ! Putative coordinates of each boundary point (initial approximation)

  real(kind=8) :: put(1:2,np)
  data put(1:2,  1) / 11.70d0,  5.70d0 /, &
       put(1:2,  2) / 12.10d0,  5.60d0 /, &
       put(1:2,  3) / 13.50d0,  8.30d0 /, &
       put(1:2,  4) / 13.10d0,  8.40d0 /, &
       put(1:2,  5) / 13.70d0,  8.60d0 /, &
       put(1:2,  6) / 13.20d0,  9.00d0 /, &
       put(1:2,  7) / 13.60d0,  9.40d0 /, &
       put(1:2,  8) / 13.40d0, 10.00d0 /, &
       put(1:2,  9) / 12.50d0, 10.00d0 /, &
       put(1:2, 10) / 13.10d0, 10.40d0 /, &
       put(1:2, 11) / 11.70d0,  6.40d0 /, &
       put(1:2, 12) / 11.60d0, 10.60d0 /, &
       put(1:2, 13) / 11.60d0, 11.20d0 /, &
       put(1:2, 14) / 12.00d0,  9.90d0 /, &
       put(1:2, 15) / 12.70d0, 11.20d0 /, &
       put(1:2, 16) / 11.80d0, 11.60d0 /, &
       put(1:2, 17) / 13.80d0, 13.50d0 /, &
       put(1:2, 18) / 13.50d0, 13.20d0 /, &
       put(1:2, 19) / 12.90d0, 14.10d0 /, &
       put(1:2, 20) / 12.80d0, 13.70d0 /, &
       put(1:2, 21) / 12.40d0, 13.30d0 /, &
       put(1:2, 22) / 11.90d0, 13.30d0 /, &
       put(1:2, 23) / 11.70d0, 12.70d0 /, &
       put(1:2, 24) / 11.60d0, 14.60d0 /, &
       put(1:2, 25) / 11.10d0, 12.80d0 /, &
       put(1:2, 26) / 10.80d0, 12.50d0 /, &
       put(1:2, 27) / 10.50d0, 12.50d0 /, &
       put(1:2, 28) / 10.70d0, 13.20d0 /, &
       put(1:2, 29) / 10.90d0, 14.20d0 /, &
       put(1:2, 30) / 10.80d0, 14.10d0 /, &
       put(1:2, 31) /  9.20d0, 15.00d0 /, &
       put(1:2, 32) /  9.40d0, 15.30d0 /, &
       put(1:2, 33) /  6.60d0, 17.60d0 /, &
       put(1:2, 34) /  8.40d0, 16.60d0 /, &
       put(1:2, 35) /  6.40d0, 20.00d0 /, &
       put(1:2, 36) / 12.50d0, 19.50d0 /, &
       put(1:2, 37) / 13.80d0, 20.60d0 /, &
       put(1:2, 38) / 11.80d0, 22.40d0 /, &
       put(1:2, 39) / 11.20d0, 20.50d0 /, &
       put(1:2, 40) / 10.00d0, 22.20d0 /, &
       put(1:2, 41) / 11.50d0, 23.00d0 /, &
       put(1:2, 42) / 11.00d0, 24.60d0 /, &
       put(1:2, 43) /  5.80d0, 23.80d0 /, &
       put(1:2, 44) /  3.50d0, 23.80d0 /, &
       put(1:2, 45) /  3.00d0, 21.50d0 /, &
       put(1:2, 46) /  5.50d0, 22.00d0 /, &
       put(1:2, 47) / 15.70d0, 12.20d0 /, &
       put(1:2, 48) / 15.30d0, 10.80d0 /, &
       put(1:2, 49) / 10.40d0, 16.20d0 /, &
       put(1:2, 50) / 11.00d0, 15.80d0 /, &
       put(1:2, 51) / 11.30d0, 16.00d0 /, &
       put(1:2, 52) / 10.40d0, 16.30d0 /, &
       put(1:2, 53) / 11.60d0, 12.00d0 /, &
       put(1:2, 54) / 13.00d0,  9.30d0 /, &
       put(1:2, 55) / 12.00d0, 13.80d0 /, &
       put(1:2, 56) /  6.00d0, 18.70d0 /, &
       put(1:2, 57) / 11.50d0,  5.70d0 /, &
       put(1:2, 58) / 11.80d0,  5.50d0 /, &
       put(1:2, 59) / 11.40d0, 10.60d0 /, &
       put(1:2, 60) / 10.70d0, 17.40d0 /, &
       put(1:2, 61) / 10.60d0, 16.60d0 /, &
       put(1:2, 62) / 10.40d0, 17.20d0 /, &
       put(1:2, 63) / 11.00d0, 19.10d0 /, &
       put(1:2, 64) / 10.00d0, 20.00d0 /, &
       put(1:2, 65) / 13.50d0,  8.80d0 /, &
       put(1:2, 66) / 13.20d0,  7.70d0 /, &
       put(1:2, 67) / 12.60d0,  7.70d0 /, &
       put(1:2, 68) / 14.20d0,  9.70d0 /, &
       put(1:2, 69) / 12.30d0,  6.70d0 /, &
       put(1:2, 70) / 13.10d0,  9.60d0 /, &
       put(1:2, 71) / 12.20d0,  5.20d0 /, &
       put(1:2, 72) / 12.50d0,  5.30d0 /, &
       put(1:2, 73) / 12.20d0,  5.40d0 /, &
       put(1:2, 74) / 12.10d0,  5.60d0 /, &
       put(1:2, 75) / 12.20d0,  5.40d0 /, &
       put(1:2, 76) / 15.00d0, 10.00d0 /, &
       put(1:2, 77) / 11.90d0,  5.30d0 /, &
       put(1:2, 78) /  2.70d0, 23.00d0 /, &
       put(1:2, 79) /  4.00d0, 22.00d0 /, &
       put(1:2, 80) /  4.50d0, 24.00d0 /, &
       put(1:2, 81) /  9.00d0, 23.50d0 /, &
       put(1:2, 82) /  7.30d0, 23.70d0 /, &
       put(1:2, 83) /  6.00d0, 21.00d0 /, &
       put(1:2, 84) / 12.00d0, 19.50d0 /, &
       put(1:2, 85) /  8.50d0, 16.00d0 /, &
       put(1:2, 86) /  7.70d0, 15.80d0 /, &
       put(1:2, 87) /  9.30d0, 15.50d0 /, &
       put(1:2, 88) /  9.90d0, 16.00d0 /, &
       put(1:2, 89) /  9.00d0, 15.60d0 /, &
       put(1:2, 90) / 10.20d0, 15.30d0 /, &
       put(1:2, 91) / 10.20d0, 14.50d0 /, &
       put(1:2, 92) /  9.70d0, 14.70d0 /, &
       put(1:2, 93) / 10.00d0, 14.40d0 /, &
       put(1:2, 94) /  9.30d0, 17.20d0 /, &
       put(1:2, 95) / 11.50d0, 14.00d0 /, &
       put(1:2, 96) / 11.30d0, 17.80d0 /, &
       put(1:2, 97) / 11.40d0, 18.40d0 /, &
       put(1:2, 98) / 10.50d0, 16.60d0 /, &
       put(1:2, 99) / 13.10d0, 19.50d0 /, &
       put(1:2,100) / 12.40d0, 19.90d0 /, &
       put(1:2,101) / 11.20d0, 20.20d0 /, &
       put(1:2,102) / 10.00d0, 21.50d0 /, &
       put(1:2,103) /  9.50d0, 16.00d0 /, &
       put(1:2,104) /  7.10d0, 16.30d0 /, &
       put(1:2,105) /  6.80d0, 17.50d0 /, &
       put(1:2,106) /  6.90d0, 17.50d0 /, &
       put(1:2,107) /  7.80d0, 15.90d0 /, &
       put(1:2,108) /  6.80d0, 17.60d0 /, &
       put(1:2,109) /  7.40d0, 17.40d0 /, &
       put(1:2,110) /  7.70d0, 17.50d0 /, &
       put(1:2,111) /  8.10d0, 17.10d0 /, &
       put(1:2,112) /  8.30d0, 17.30d0 /, &
       put(1:2,113) /  8.50d0, 15.30d0 /, &
       put(1:2,114) /  8.90d0, 15.30d0 /, &
       put(1:2,115) /  9.70d0, 15.60d0 /, &
       put(1:2,116) /  7.20d0, 16.30d0 /, &
       put(1:2,117) /  7.40d0, 16.60d0 /, &
       put(1:2,118) /  6.80d0, 16.80d0 /, &
       put(1:2,119) /  7.00d0, 16.80d0 /, &
       put(1:2,120) /  6.10d0, 19.50d0 /, &
       put(1:2,121) /  6.20d0, 18.20d0 /, &
       put(1:2,122) / 12.50d0, 14.50d0 /, &
       put(1:2,123) / 12.10d0, 14.50d0 /, &
       put(1:2,124) /  2.75d0, 22.15d0 /, &
       put(1:2,125) / 14.60d0, 12.80d0 /, &
       put(1:2,126) / 13.00d0, 13.20d0 /, &
       put(1:2,127) / 10.20d0, 23.30d0 /, &
       put(1:2,128) / 10.70d0, 11.25d0 /, &
       put(1:2,129) / 15.30d0, 11.30d0 /, &
       put(1:2,130) / 11.40d0,  7.50d0 /, &
       put(1:2,131) / 11.65d0,  9.00d0 /, &
       put(1:2,132) / 11.60d0,  7.25d0 /

  ! 'Center' of each country

  real(kind=8) :: w(2,nstates)

  ! Target values

  ! Population taken from
  ! http://en.wikipedia.org/wiki/List_of_countries_by_population,
  ! accessed on May 21th, 2013.

  ! Area in km^2 taken from
  ! https://en.wikipedia.org/wiki/List_of_countries_and_dependencies_by_area
  ! accessed on May 21th, 2013.

  ! GDP in millions of $US taken from:
  ! http://en.wikipedia.org/wiki/List_of_countries_by_GDP_%28nominal%29
  ! accessed on May 21th, 2013.

  ! Particular cases:
  !
  ! Alaska population and area: https://en.wikipedia.org/wiki/Alaska, 
  ! accessed on May 21th, 2013.

  ! Alaska GDP from: http://en.wikipedia.org/wiki/List_of_U.S._states_by_GDP,
  ! accessed on May 21th, 2013.

  ! French Guiana population and GDP: https://en.wikipedia.org/wiki/French_Guiana, 
  ! accessed on May 21th, 2013.

  ! Population:
  !
  !  1 Cuba            11,163,934
  !  2 Canada          35,056,064
  !  3 USA            315,150,551 (USA 315,882,000 - Alaska 731,449)
  !  4 Mexico         112,336,538
  !  5 America Central 44,463,381 (Guatemala 15,438,384, Hunduras 8,385,072, Nicaragua 6,071,045, El Salvador 6,183,000, 
  !                                Costa Rica 4,667,096, Panama 3,405,813, and Belize 312,971)
  !  6 Colombia        47,059,000
  !  7 Venezuela       28,946,101
  !  8 Guianas          1,548,123 (Guyana 784,894, Suriname 534,189, and French Guiana 229,040)
  !  9 Brasil         193,946,886
  ! 10 Ecuador         15,489,300
  ! 11 Peru            30,475,144
  ! 12 Bolivia         10,389,913
  ! 13 Paraguay         6,672,631
  ! 14 Chile           16,634,603
  ! 15 Argentina       40,117,096
  ! 16 Alaska             731,449
  ! 17 Uruguay          3,286,314

!!$  real(kind=8) :: tarval(nstates)
!!$  data tarval(1:nstates) / 11163934.0d0, &
!!$                           35056064.0d0, &
!!$                          315150551.0d0, &
!!$                          112336538.0d0, &
!!$                           44463381.0d0, &
!!$                           47059000.0d0, &
!!$                           28946101.0d0, &
!!$                            1548123.0d0, &
!!$                          193946886.0d0, &
!!$                           15489300.0d0, &
!!$                           30475144.0d0, &
!!$                           10389913.0d0, &
!!$                            6672631.0d0, &
!!$                           16634603.0d0, &
!!$                           40117096.0d0, &
!!$                             731449.0d0, &
!!$                            3286314.0d0 /

  ! Areas:
  !
  !  1 Cuba            109,884
  !  2 Canada        9,984,670
  !  3 USA           7,911,237 ( USA 9,629,091 - Alaska 1,717,854 )
  !  4 Mexico        1,964,375
  !  5 America Central 476,278 (Guatemala 108,889, Hunduras 112,492, Nicaragua 130,373, El Salvador 21,041, 
  !                             Costa Rica 51,100, Panama 75,417, and Belize 22,966) 
  !  6 Colombia      1,141,748
  !  7 Venezuela       912,050
  !  8 Guianas         462,323 (Guyana 214,969, Suriname 163,820, and French Guiana 83,534)
  !  9 Brasil        8,514,877
  ! 10 Ecuador         256,369
  ! 11 Peru          1,285,216
  ! 12 Bolivia       1,098,581
  ! 13 Paraguay        406,752
  ! 14 Chile           756,102
  ! 15 Argentina     2,780,400
  ! 16 Alaska        1,717,854
  ! 17 Uruguay         181,034

  real(kind=8) :: tarval(nstates)
  data tarval(1:nstates) /  109884.0d0, &
                           9984670.0d0, &
                           7911237.0d0, &
                           1964375.0d0, &
                            476278.0d0, &
                           1141748.0d0, &
                            912050.0d0, &
                            462323.0d0, &
                           8514877.0d0, &
                            256369.0d0, &
                           1285216.0d0, &
                           1098581.0d0, &
                            406752.0d0, &
                            756102.0d0, &
                           2780400.0d0, &
                           1717854.0d0, &
                            181034.0d0 /

!!$  ! Original areas in a previous version of this problem (unkonwn source)
!!$
!!$  real(kind=8) :: tarval(nstates)
!!$  data tarval(1:nstates) / 115.0d0, 9922.0d0, 7885.0d0, 1973.0d0,  543.0d0, 1139.0d0, &
!!$                           912.0d0,  470.0d0, 8512.0d0,  461.0d0, 1285.0d0, 1099.0d0, &
!!$                           407.0d0,  752.0d0, 2767.0d0, 1478.0d0,  187.0d0 /

  ! GDP (gross domestic product) in millions of $US:
  !
  !  1 Cuba                68,715
  !  2 Canada           1,736,869
  !  3 USA             14,991,300
  !  4 Mexico           1,155,206
  !  5 America Central    167,854 (Guatemala 46,898, Hunduras 17,447, Nicaragua 7,297, El Salvador 23,054, 
  !                                Costa Rica 41,007, Panama 30,677, and Belize 1,474)
  !  6 Colombia           333,185
  !  7 Venezuela          315,893
  !  8 Guianas             12,126 (Guyana 2,577, Suriname 4,610, and French Guiana 4,939)
  !  9 Brasil           2,476,651
  ! 10 Ecuador             66,381
  ! 11 Peru               180,464
  ! 12 Bolivia             23,949
  ! 13 Paraguay            22,890
  ! 14 Chile              248,592
  ! 15 Argentina          485,416
  ! 16 Alaska              45,600
  ! 17 Uruguay             46,710

!!$  real(kind=8) :: tarval(nstates)
!!$  data tarval(1:nstates) /     68715.0d0, &
!!$                             1736869.0d0, &
!!$                            14991300.0d0, &
!!$                             1155206.0d0, &
!!$                              167854.0d0, &
!!$                              333185.0d0, &
!!$                              315893.0d0, &
!!$                               12126.0d0, &
!!$                             2476651.0d0, &
!!$                               66381.0d0, &
!!$                              180464.0d0, &
!!$                               23949.0d0, &
!!$                               22890.0d0, &
!!$                              248592.0d0, &
!!$                              485416.0d0, &
!!$                               45600.0d0, &
!!$                               46710.0d0 /

  ! Colors for the graphical representation

  character * 6 color(nstates)
  data color( 1) /'yellow'/, &
       color( 2) /'yellow'/, &
       color( 3) /'purple'/, &
       color( 4) /'green '/, &
       color( 5) /'coral '/, &
       color( 6) /'orange'/, &
       color( 7) /'green '/, &
       color( 8) /'purple'/, &
       color( 9) /'yellow'/, &
       color(10) /'yellow'/, &
       color(11) /'coral '/, &
       color(12) /'orange'/, &
       color(13) /'purple'/, &
       color(14) /'yellow'/, &
       color(15) /'green '/, &
       color(16) /'purple'/, &
       color(17) /'orange'/

end module modamerica
