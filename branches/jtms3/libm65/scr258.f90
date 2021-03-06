subroutine scr258(isync,idat,ndir,ichan)

  integer*1 isync(43)
  integer*1 idat(215)
  integer*1 ichan(258)

  integer indx(258)
  data indx/                                            &
      -1,   1, 129,  65, 193,  33,  -2, 161,  97,  17,   &  ! 10
     145,  81,  -3, 209,  49, 177, 113,   9,  -4, 137,   &  ! 20
      73, 201,  41, 169,  -5, 105,  25, 153,  89,  57,   &  ! 30
      -6, 185, 121,   5, 133,  69,  -7, 197,  37, 165,   &  ! 40
     101,  21,  -8, 149,  85, 213,  53, 181,  -9, 117,   &  ! 50
      13, 141,  77, 205, -10,  45, 173, 109,  29, 157,   &  ! 60
     -11,  93,  61, 189, 125,   3, -12, 131,  67, 195,   &  ! 70
      35, 163, -13,  99,  19, 147,  83, 211, -14,  51,   &  ! 80
     179, 115,  11, 139, -15,  75, 203,  43, 171, 107,   &  ! 90
     -16,  27, 155,  91,  59, 187, -17, 123,   7, 135,   &  !100
      71, 199, -18,  39, 167, 103,  23, 151, -19,  87,   &  !110
     215,  55, 183, 119, -20,  15, 143,  79, 207,  47,   &  !120
     -21, 175, 111,  31, 159,  95, -22,  63, 191, 127,   &  !130
       2, 130, -23,  66, 194,  34, 162,  98, -24,  18,   &  !140
     146,  82, 210,  50, -25, 178, 114,  10, 138,  74,   &  !150
     -26, 202,  42, 170, 106,  26, -27, 154,  90,  58,   &  !160
     186, 122, -28,   6, 134,  70, 198,  38, -29, 166,   &  !170
     102,  22, 150,  86, -30, 214,  54, 182, 118,  14,   &  !180
     -31, 142,  78, 206,  46, 174, -32, 110,  30, 158,   &  !190
      94,  62, -33, 190, 126,   4, 132,  68, -34, 196,   &  !200
      36, 164, 100,  20, -35, 148,  84, 212,  52, 180,   &  !210
     -36, 116,  12, 140,  76, 204, -37,  44, 172, 108,   &  !220
      28, 156, -38,  92,  60, 188, 124,   8, -39, 136,   &  !230
      72, 200,  40, 168, -40, 104,  24, 152,  88,  56,   &  !240
     -41, 184, 120,  16, 144,  80, -42, 208,  48, 176,   &  !250
     112,  32, -43, 160,  96,  64, 192, 128/
  save

  if(ndir.gt.0) then
     do i=1,258
        j=indx(i)
        if(j.lt.0) ichan(i)=isync(-j)
        if(j.gt.0)  ichan(i)=idat(j)
     enddo
  else
     do i=1,258
        j=indx(i)
        if(j.lt.0) isync(-j)=ichan(i)
        if(j.gt.0) idat(j)=ichan(i)
     enddo
  endif

end subroutine scr258
