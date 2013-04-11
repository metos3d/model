
#ifndef DIC_OPTIONS_H
#define DIC_OPTIONS_H

#define _RL Real*8
#define _RS Real*8
#define _d D

#define ALLOW_PTRACERS
#define DIC_BIOTIC
#define DIC_NO_NEG

#define _hFacC(i,j,k,bi,bj) hFacC(i,j,k,bi,bj)
#define _recip_hFacC(i,j,k,bi,bj) recip_hFacC(i,j,k,bi,bj)

#endif