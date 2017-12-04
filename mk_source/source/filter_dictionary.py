import numpy as np

filters = {}

filters[2190] ={ 'name'  : 'K',
                 'lambda': 2190e-9,
                 'active': 1,
                 'filename': ['filter_measures_K.dat']
             }

filters[445]  ={ 'name'  : 'B',
                 'lambda': 445e-9,
                 'active': 1,
                 'filename': ['filter_measures_B.dat']
             }

filters[365] ={ 'name'  : 'U',
                 'lambda': 365e-9,
                 'active': 1,
                 'filename': ['filter_measures_U.dat']
             }

filters[337] ={ 'name'  : 'U_HST',
                 'lambda': 337.5e-9,
                 'active': 0,
                 'filename': ['filter_measures_F336W.dat']
             }

filters[1005] ={ 'name'  : 'y',
                 'lambda': 1005e-9,
                 'active': 0,
                 'filename': ['filter_measures_y.dat']
             }

filters[1020] ={ 'name'  : 'Y',
                 'lambda': 1020e-9,
                 'active': 0,
                 'filename': ['filter_measures_Y.dat']
             }

filters[354]  ={ 'name'  : 'u',
                 'lambda': 354e-9,
                 'active': 0,
                 'filename': ['filter_measures_u.dat']
             }

filters[475]  ={ 'name'  : 'g',
                 'lambda': 475e-9,
                 'active': 1,
                 'filename': ['filter_measures_g.dat','filter_measures_F475W.dat']
             }

#filters[520]  ={ 'name'  : 'g',
#                 'lambda': 529e-9,
#                 'active': 1,
#                 'filename': ['filter_measures_g.dat','filter_measures_F475W.dat']
#             }

filters[850]  ={ 'name'  : 'z',
                 'lambda': 850e-9,
                 'active': 1,
                 'filename': ['filter_measures_z.dat','filter_measures_F850W.dat']
             }

filters[551]  ={ 'name'  : 'V',
                 'lambda': 551e-9,
                 'active': 1,
                 'filename': ['filter_measures_V.dat']
             }

filters[596]  ={ 'name'  : 'V_HST',
                 'lambda': 595.6e-9,
                 'active': 0,
                 'filename': ['filter_measures_F606W.dat']
             }


filters[1220] ={ 'name'  : 'J',
                 'lambda': 1220e-9,
                 'active': 1,
                 'filename': ['filter_measures_J.dat']
             }

filters[10469] ={ 'name'  : 'J1',
                 'lambda': 1046.9e-9,
                 'active': 0,
                 'filename': ['filter_measures_J1.dat']
             }

filters[1150] ={ 'name'  : 'J_HST',
                 'lambda': 1150e-9,
                 'active': 0,
                 'filename': ['filter_measures_F110W.dat']
             }

filters[2150] ={ 'name'  : 'Ks',
                 'lambda': 2150e-9,
                 'active': 1,
                 'filename': ['filter_measures_Ks.dat']
             }

filters[1630] ={ 'name'  : 'H',
                 'lambda': 1630e-9,
                 'active': 1,
                 'filename': ['filter_measures_H.dat']
             }

filters[1545] ={ 'name'  : 'H_HST',
                 'lambda': 1545e-9,
                 'active': 0,
                 'filename': ['filter_measures_F160W.dat']
             }

filters[658]  ={ 'name'  : 'R',
                 'lambda': 658e-9,
                 'active': 1,
                 'filename': ['filter_measures_R.dat']
             }

filters[625]  ={ 'name'  : 'r',
                 'lambda': 625e-9,
                 'active': 1,
                 'filename': ['filter_measures_r.dat','filter_measures_F625W.dat']
             }

filters[12130]={ 'name'  : 'C',
                 'lambda': 12130e-9,
                 'active': 0,
                 'filename': ['filter_measures_C.dat']
             }

filters[806]  ={ 'name'  : 'I',
                 'lambda': 806e-9,
                 'active': 1,
                 'filename': ['filter_measures_I.dat','filter_measures_F814W.dat']
               }

filters[775]  ={ 'name'  : 'i',
                 'lambda': 775e-9,
                 'active': 1,
                 'filename': ['filter_measures_i.dat','filter_measures_F775W.dat']
               }
