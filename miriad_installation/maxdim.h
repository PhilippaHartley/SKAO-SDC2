C=======================================================================
C $Id: maxdim.h,v 1.15 2020/05/18 05:33:42 mirmgr Exp $
C-----------------------------------------------------------------------
C     Size of an INTEGER array used to implement a memory heap.  This
C     array is the sole variable in blank COMMON in Miriad.  Trial-and-
C     error compilations on an x86-64 system with gcc/g77 show that the
C     limit on MAXBUF for which most tasks build successfully is
C     1073741823 (2**30 - 1) which corresponds to 4GiB.  With MAXDIM
C     less than 32768 the limit for which all tasks build successfully
C     is about 260000000; unsuccessful links produce messages about
C     truncated relocations, imom being the worst offender.
C     The default value allocates 128MiB (for normal 4-byte INTEGERs).
      INTEGER   MAXBUF
C      PARAMETER(MAXBUF = 1073741823)
      PARAMETER(MAXBUF = 32*1024*1024)

C     Maximum image axis length.  Array dimensions are typically a few
C     times MAXDIM (never MAXDIM**2) so MAXDIM is associated with a much
C     smaller allocation of static memory than MAXBUF.  Thus the default
C     value of MAXDIM is quite generous. 
C     The current value of MAXDIM is the largest value that keeps the
C     memory allocation to store a single image plane below the current
C     limit of 8GB. Note that invert is limited to 32768, the highest
C     power of 2 below the limit.
      INTEGER   MAXDIM
      PARAMETER(MAXDIM = 65536)

C     Maximum number of antennas (ATA=64).
      INTEGER   MAXANT
      PARAMETER(MAXANT = 1024)

C     Maximum number of baselines, including autocorrelations.
      INTEGER   MAXBASE
      PARAMETER(MAXBASE = ((MAXANT*(MAXANT+1))/2))

C     Maximum number of channels in spectral data.
      INTEGER   MAXCHAN
      PARAMETER(MAXCHAN = 70000)

C     Maximum number of windows in visibility data.
      INTEGER   MAXWIN
      PARAMETER(MAXWIN = 48)

C     Maximum number of wideband channels.
      INTEGER   MAXWIDE
      PARAMETER(MAXWIDE = 18)

C     Maximum number of frequency bins in gains table
      INTEGER   MAXFBIN
      PARAMETER(MAXFBIN = 16)

C     Maximum number of mosaic pointings
      INTEGER   MAXPNT
      PARAMETER(MAXPNT = 20000)
C=======================================================================