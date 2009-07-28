C
      PROGRAM CKINTP
C
C----------------------------------------------------------------------C
C     VERSION 3.6
C     CHANGES FROM VERSION 1.0
C     1.  Changed from REAL*8 to DOUBLE PRECISION
C     CHANGES FROM VERSION 1.1
C     1.  Changed CHARACTER*100 to CHARACTER*80
C     2.  Added THERMO "ALL" option
C     3.  Write LENICK, LENRCK, LENCCK to binary file
C     4.  Allow reaction species to end in '=' or '-'
C     5.  Allow real values of elemental composition in THERMO cards
C     6.  Allow upper/lower case input
C     CHANGES FROM VERSION 1.2
C     1.  Reaction delimiters are now "=" or "<=>" if reversible,
C                                            " =>" if irreversible.
C     2.  Fixed an error with IFIRCH(LINE) in IPPLEN
C     CHANGES FROM VERSION 1.3
C     1.  Add "unix" change blocks
C     CHANGES FROM VERSION 1.4
C     1.  Modify OPEN statements
C     CHANGES FROM VERSION 1.5
C     1.  Correct molecules to moles unit conversion
C     2.  Correct UPCASE to avoid dimensioning errors
C     CHANGES FROM VERSION 1.7
C     1.  Further correction of molecules conversion for fall-off
C         and third-body reactions
C     CHANGES FOR VERSION 1.8
C     1.  Change Subroutine CKUNIT to parse LINE instead of SUB(*)
C         in order to correct misinterpretation of unit strings
C         with slashes.
C     CHANGES FOR VERSION 1.9
C     1.  First record of binary file now consists of a character
C         string version, precision, and logical error flag
C     CHANGES FOR VERSION 2.0
C     1.  Error in UPCASE could cause interpreter to ignore some
C         keywords.
C     CHANGES FOR VERSION 2.1
C     1.  10/18/90 (F. Rupley):
C         Error in scaling the pre-exponential constants RPAR(3,*)
C         where REV is declared, and FPAL(3,*) for fall-off reactions,
C         as RPAR(3,II)*EFAC should read RPAR(3,NREV), and
C            FPAL(3,II)*EFAC should read FPAL(3,NFAL).
C         This error was introduced in CKINTERP.15 during refinement
C         Dof units conversion routines.
C     2.  Subroutine CKDUP modified to recognize that two reactions
C         may be duplicate except for a third-body species in a
C         fall-off reaction.
C     CHANGES FOR VERSION 2.2
C     1.  11/14/90 (F. Rupley per M. Coltrin):
C         Initialize variable NCHRG
C      CHANGES FOR VERSION 2.3
C     1.  In CKPREAC, error correction of 10/18/90 (above, V2.1).
C     CHANGES FOR VERSION 2.4
C     1.  Additional checking of TLO,TMID,THI for species -
C         a) set initial values at -1.
C         b) if user has not provided a TLO,TMID, or THI, use the
C            values provided by THERMO.DAT.
C         c) check that TLO < THI, TLO <= TMID <= THI
C     CHANGES FOR VERSION 2.5
C     1.  Need to get TLO,THI,TMID from database BEFORE reading
C         user's THERMO data (unless THERMO ALL option is used)
C     CHANGES FOR VERSION 2.6
C     1.  LENRCK lengthened by II+NREV to reflect additional
C         work space needed by CKRAT for a 4th parameter
C         (perturbation factor).
C     CHANGES FOR VERSION 2.7
C     1.  Two otherwise duplicate reactions are unique if one
C         is a third body reaction and the other not.
C     CHANGES FOR VERSION 2.8
C     1.  Change output format to print all 16 characters for
C         a species name.
C     CHANGES FOR VERSION 2.9 (2/24/92 F. Rupley)
C     1.  Check that reverse (REV) parameters were given when
C         RTL reverse Teller-Landauer parameters are given.
C     2.  Add 2*II to length of real work space
C     CHANGES FOR VERSION 3.0 (4/13/92 F. Rupley per M. Coltrin)
C     1.  Correct logic in CKDUP, add argument to call list.
C     CHANGES FOR VERSION 3.1 (2/24/93 F. Rupley per C. Westbrook,LLNL)
C     1.  Problem in CKREAC for species starting with "M", where
C         "+M" is signal for third-body.
C     CHANGES FOR VERSION 3.2 (11/11/93 F. Rupley per T.U.Delft)
C     1.  Ensure that SUBROUTINE CKUNIT does not check for units beyond
C         end of LINE.
C     CHANGES FOR VERSION 3.3 (1/26/94 F. Rupley per R. Kee)
C     1.  Real stoichometric coefficients used in a supplemental way;
C         NRNU total number of reactions using real stoichometry,
C         IRNU array of reaction numbers, RNU real coefficients.
C     CHANGES FOR VERSION 3.4 (3/15/94 F. Rupley)
C     1.  DOS/PC compatibility effort includes adding file names to
C         OPEN statements, removing unused variables in CALL lists,
C         unusued but possibly initialized variables.
C     CHANGES FOR VERSION 3.5 (4/19/94 F. Rupley)
C     1.  Fix bug with index KSPEC(N) for CKBAL and CKRBAL.
C
C     CKINTP interprets a formatted ASCII representation of a
C     chemical reaction mechanism and creates the binary file LINK
C     required by CHEMKIN.  CKINTP is dimensioned as follows:
C
C     MDIM = maximum number of elements in a problem;             (10)
C     KDIM = maximum number of species in a problem;             (100)
C     MAXTP= maximum number of temperatures used to fit            (3)
C            thermodynamic properties of species
C     NPC  = number of polynomial coefficients to fits             (5)
C     NPCP2= number of fit coefficients for a temperature range    (7)
C     IDIM = maximum number of reactions in a mechanism;         (500)
C     NPAR = number of Arrhenius parameters in a reaction;         (3)
C     NLAR = number of Landau-Teller parameters in a reaction;     (2)
C     NFAR = number of fall-off parameters in a reaction;          (8)
C     MAXSP= maximum number of species in a reaction               (6)
C     MAXTB= maximum number of third bodies for a reaction        (10)
C     LSYM = character string length of element and species names (16)
C
C     User input is read from LIN (Unit15), a thermodynamic database
C     is read from LTHRM (Unit17), printed output is assigned to LOUT
C     (Unit16), and binary data is written to LINC (Unit25).
C
C     REQUIRED ELEMENT INPUT: (Subroutine CKCHAR)          (DIMENSION)
C
C        The word 'ELEMENTS' followed by a list of element
C        names, terminated by the word 'END';
C
C        The resulting element data stored in LINK is:
C        MM       - integer number of elements found
C        ENAME(*) - CHARACTER*(*) array of element names        (MDIM)
C        AWT(*)   - real array of atomic weights;               (MDIM)
C                   default atomic weights are those on
C                   atomic weight charts; if an element
C                   is not on the periodic chart, or if
C                   it is desirable to alter its atomic
C                   weight, this value must be included
C                   after the element name, enclosed by
C                   slashed, i.e., D/2.014/
C
C     REQUIRED SPECIES INPUT: (Subroutine CKCHAR)
C
C        The word 'SPECIES' followed by a list of species
C        names, terminated by the word 'END';
C
C        The resulting species data stored in LINK is:
C        KK       - integer number of species found
C        KNAME(*) - CHARACTER*(*) array of species names        (KDIM)
C
C     OPTIONAL THERMODYNAMIC DATA: (Subroutine CKTHRM)
C     (If this feature is not used, thermodynamic properties are
C     obtained from a CHEMKIN database.)  The format for this option
C     is the word 'THERMO' followed by any number of 4-line data sets:
C
C     Line 1: species name, optional comments, elemental composition,
C             phase, T(low), T(high), T(mid), additional elemental
C             composition, card number (col. 80);
C             format(A10,A14,4(A2,I3),A1,E10.0,E10.0,E8.0,(A2,I3),I1)
C     Line 2: coefficients a(1--5) for upper temperature range,
C             card number (col. 80);
C             format(5(e15.0),I1)
C     Line 3: coefficients a(6--7) for upper temperature range,
C             coefficients a(1--3) for lower temperature range,
C             card number (col. 80);
C             format(5(e15.0),I1)
C     Line 4: coefficients a(4--7) for lower temperature range,
C             card number (col. 80);
C             format(4(e15.0),I1)
C
C     End of THERMO data is indicated by 'END' line or new keyword.
C
C        The resulting thermodynamic data stored in LINK are:
C        WTM(*)   - real array of molecular weights             (KDIM)
C        KNCF(*,*)- integer composition of species         (MDIM,KDIM)
C        KPHSE(*) - integer phase of a species;                 (KDIM)
C                   -1(solid), 0(gas), +1(liquid).
C        KCHRG(*) - ionic charge of a species;                  (KDIM)
C                   = 0 except in presence/absence of electrons
C                   = +n in absence of n electrons
C                   = -n in presence of n electons
C        NCHRG    - integer number of species with KCHRG<>0
C        NT(*)    - array of number of temperatures used        (KDIM)
C                   in fits
C        T(*,*)   - array of temperatures used in fits    (MAXTP,KDIM)
C        A(N,L,K) - Thermodynamic properties for      (NPC+2,NTR,KDIM)
C                   species K consists of polynomial
C                   coefficients for fits to
C                   CP/R = SUM (A(N,L,K)*Temperature**(N-1), N=1,NPC+2)
C                          where  T(L,K) <= Temperature < T(L+1,K),
C                   and,
C                   N=NPC+1 is formation enthalpy HO/R = A(NPC+1,L,K),
C                   N=NPC+2 is formation entropy  SO/R = A(NPC+2,L,K)
C
C     OPTIONAL REACTION INPUT:
C     Reaction data is input after all ELEMENT, SPECIES and THERMO
C     data in the following format:
C
C     1) (Subroutine CKREAC)
C        The first line contains the keyword 'REACTIONS' and an
C        optional description of units:
C
C           'MOLES' - (default), pre-exponential units are moles-sec-K;
C           'MOLECULES' - pre-exponential units are molecules and
C                         will be converted to moles.
C           'KELVINS' - activation energies are Kelvins, else the
C                       activation energies are converted to Kelvins;
C           'CAL/MOLE' - (default), activation energies are cal/mole;
C           'KCAL/MOLE' - activation energies are Kcal/mole;
C           'JOULES/MOLE' - activation energies are joules/mole;
C           'KJOULES/MOLE' - activation energies are Kjoules/mole.
C
C        A description of each reaction is expected to follow.
C        Required format for a reaction is a list of '+'-delimited
C        reactants, followed by a list of '+'-delimited reactants,
C        each preceded by its stoichiometric coefficient if greater
C        than 1;  separating the reactants from the products is a '='
C        if reversible reaction, else a '=>'.  Following the reaction
C        string on the same line are the space-delimited Arrhenius
C        coefficients.
C
C        If the reaction contains a third body, this is indicated by
C        by the presence of an 'M' as a reactant or product or both,
C        and enhancement factors for third-bodies may be defined on
C        additional lines as described in (2).
C
C        If the reaction contains a radiation wavelength, this is
C        indicated by the presence of an 'HV' either as a reactant
C        or as a product.  Unless otherwise defined on additional
C        lines as described in (2), the value of the wavelength is
C        -1.0 if a reactant or +1.0 if a product.
C
C        If the reaction is a fall-off reaction, this is indicated
C        either by a '(+M)' or a '(+KNAME(K))', and there must be
C        additional lines as described in (2) to define fall-off
C        parameters.
C
C    2)  (Subroutine CKAUXL)
C        Additional information for a reaction is given on lines
C        immediately following the reaction description; this data
C        will consist of a 'keyword' to denote the type of data,
C        followed by a '/', then the required parameters for the
C        keyword, followed by another '/'.  There may be more than
C        one keyword per line, and there may be any number of lines.
C        The keywords and required parameters are as follows:
C
C        KNAME(K)/efficiency value/ - species (K) is an enhanced
C                third body in the reaction
C        HV/wavelength/ - radiation wavelength parameter
C        LT/val1 val2/ - Landau-Teller coefficients
C        LOW/val1 val2 val3/ - low fall-off parameters
C        TROE/val1 val2 val3 val4/ - Troe fall-off parameters;
C                                    if val4 is omitted, a default
C                                    parameter will be used
C        SRI/val1 val2 val3 val4/ - SRI fall-off parameters;
C                                   if val4 is omitted, a default
C                                   parameter will be used
C           (it is an error to have both LT and Fall-off defined)
C        REV/par1 par2 par3/ - reverse parameters given
C        RLT/val1 val2/ - Landau-Teller coefficients for reverse
C           (it is an error if REV given and not RLT)
C
C     The end of all reaction data is indicated by an 'END' card or
C     <eof>.
C
C     Resulting reaction data stored in LINC are:
C       II        - integer number of reactions found
C       PAR(*,*)  - array of real Arrhenius coefficients   (NPAR,IDIM)
C       NSPEC(*)  - total number of species in a reaction       (IDIM)
C                   if NSPEC < 0, reaction is irreversible
C       NREAC(*)  - number of reactants only                    (IDIM)
C       NUNK(*,*) - array of species numbers for reaction (MAXSP,IDIM)
C       NU(*,*)   - array of stoichiometric coefficients  (MAXSP,IDIM)
C                   of species in a reaction, negative=reactant,
C                   positive=product
C
C       NWL       - number of reactions with radiation wavelength
C       IWL(*)    - integer reaction numbers                    (IDIM)
C       WL(*)     - real radiation wavelengths                  (IDIM)
C
C       NTHB      - number of reactions with third bodies
C       ITHB      - integer reaction numbers                    (IDIM)
C       NTBS(*)   - total number of enhanced species for NTHB   (IDIM)
C       NKTB(*,*) - species numbers of enhanced species   (MAXTB,IDIM)
C       AIK(*,*)  - enhancement factors                   (MAXTB,IDIM)
C
C       NFAL      - number of fall-off reactions
C       IFAL(*)   - integer reaction numbers                    (IDIM)
C       KFAL(*)   - integer species number for which
C                   concentrations are a factor in fall-off
C                   calculation
C       IFOP(*)   - integer fall-off type number                (IDIM)
C                   = 0 if fall-off reaction is found
C                   = 1 for Lindemann form
C                   = 2 for 6-parameter Troe form
C                   = 3 for 7-parameter Troe form
C                   = 4 for SRI form
C       PFAL(*,*) - fall-off parameters                    (NFAR,IDIM)
C
C       NLAN      - number of reactions with Landau-Teller
C       ILAN      - integer reaction numbers                    (IDIM)
C       PLAN      - Landau-Teller parameters               (NLAR,IDIM)
C
C       NREV      - number of reactions with reverse parameters
C       IREV(*)   - integer reaction numbers                    (IDIM)
C       RPAR(*,*) - parameters                             (NPAR,IDIM)
C
C       NRLT      - number of reactions with reverse parameters
C                   and Landau-Teller parameters
C       IRLT(*)   - integer reaction numbers                    (IDIM)
C       RLAN(*,*) - reverse Teller-Laudauer parameters     (NLAR,IDIM)
C
C----------------------------------------------------------------------C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      PARAMETER (MDIM=50, KDIM=500, MKDIM=MDIM*KDIM, IDIM=500, LSYM=16,
     1           NPAR=3, NPIDIM=IDIM*NPAR, NPC=5, NPCP2=NPC+2, MAXTP=3,
     2           NTR=MAXTP-1, NKTDIM=NTR*NPCP2*KDIM, MAXSP=6, MAXTB=10,
     3           NLAR=2, NSIDIM=MAXSP*IDIM, NTIDIM=MAXTB*IDIM,
     4           NLIDIM=NLAR*IDIM, NFAR=8, NFIDIM=NFAR*IDIM,
     5           NTDIM=KDIM*MAXTP, NIDIM=11*IDIM, LIN=15, LOUT=16,
     6           LTHRM=17, LINC=25, CKMIN=1.0E-3, MAXORD=10, 
     7           NOIDIM=MAXORD*IDIM)
C
      CHARACTER KNAME(KDIM)*(LSYM), ENAME(MDIM)*(LSYM), SUB(80)*80,
     1          KEY(5)*4, LINE*80, IUNITS*80, AUNITS*4, EUNITS*4,
     2          UPCASE*4, VERS*(LSYM), PREC*(LSYM)
C
      DIMENSION AWT(MDIM), KNCF(MDIM,KDIM), WTM(KDIM), KPHSE(KDIM),
     1          KCHRG(KDIM), A(NPCP2,NTR,KDIM), T(MAXTP,KDIM), NT(KDIM),
     2          NSPEC(IDIM), NREAC(IDIM), NU(MAXSP,IDIM),
     3          NUNK(MAXSP,IDIM), PAR(NPAR,IDIM), IDUP(IDIM),IREV(IDIM),
     4          RPAR(NPAR,IDIM), ILAN(IDIM), PLAN(NLAR,IDIM),
     5          IRLT(IDIM), RLAN(NLAR,IDIM), IWL(IDIM),  WL(IDIM),
     6          IFAL(IDIM), IFOP(IDIM), KFAL(IDIM), PFAL(NFAR,IDIM),
     7          ITHB(IDIM),NTBS(IDIM),AIK(MAXTB,IDIM),NKTB(MAXTB,IDIM),
     8          IRNU(IDIM), RNU(MAXSP,IDIM), IORD(IDIM), 
     9          KORD(MAXORD,IDIM), RORD(MAXORD,IDIM)
      DIMENSION VALUE(5)
C
      LOGICAL KERR, THERMO, ITHRM(KDIM)
C
C     Initialize variables
C
      DATA KEY/'ELEM','SPEC','THER','REAC','END'/, KERR/.FALSE./,
     1     ITASK,NCHRG,MM,KK,II,NLAN,NFAL,NTHB,NREV,NRLT,NWL,
     *     NRNU,NORD/13*0/,
     2     ENAME,AWT/MDIM*' ',MDIM*0.0/, THERMO/.TRUE./,
     3     T/NTDIM*-1.0/, KNAME,WTM,NT,KPHSE,KCHRG,ITHRM
     4     /KDIM*' ', KDIM*0.0, KDIM*3, KDIM*0, KDIM*0, KDIM*.FALSE./,
     5     WL,IFOP,NTBS,IDUP /IDIM*0.0, IDIM*-1, IDIM*0, IDIM*0/,
     6     NSPEC,NREAC,IREV,ILAN,IRLT,IWL,IFAL,KFAL,ITHB,IRNU,IORD
     7     /NIDIM*0/
C
      DATA NUNK,NU/NSIDIM*0, NSIDIM*0/, NKTB,AIK/NTIDIM*0,NTIDIM*-1.0/
      DATA RNU/NSIDIM*0.0/, KORD/NOIDIM*0/, RORD/NOIDIM*0.0/
      DATA PAR,RPAR/NPIDIM*0.0, NPIDIM*0.0/
      DATA PLAN,RLAN/NLIDIM*0.0, NLIDIM*0.0/
      DATA PFAL/NFIDIM*0.0/, KNCF/MKDIM*0.0/, A/NKTDIM*0.0/
C----------------------------------------------------------------------C
C
      OPEN (LOUT, FORM='FORMATTED', STATUS='UNKNOWN', FILE='chem.out')
C
      VERS = '3.6'
      WRITE  (LOUT, 15) VERS(:3)
   15 FORMAT (/
     1' CHEMKIN INTERPRETER OUTPUT: CHEMKIN-II Version ',A,' Apr. 1994'
C*****precision > double
     2/'                              DOUBLE PRECISION'/)
      PREC = 'DOUBLE'
C*****END precision > double
C*****precision > single
C     2/'                              SINGLE PRECISION'/)
C      PREC = 'SINGLE'
C*****END precision > single
C
C        START OF MECHANISM INTERPRETATION
C
      OPEN (LIN, FORM='FORMATTED', STATUS='UNKNOWN', FILE='chem.inp')
C
  100 CONTINUE
      LINE = ' '
      READ (LIN,'(A)',END=5000) LINE
  105 CONTINUE
      ILEN = IPPLEN(LINE)
      IF (ILEN .EQ. 0) GO TO 100
C
      CALL CKISUB (LINE(:ILEN), SUB, NSUB)
C
C        IS THERE A KEYWORD?
C
      CALL CKCOMP ( UPCASE(SUB(1), 4) , KEY, 5, NKEY)
      IF (NKEY .GT. 0) ITASK = 0
C
      IF (NKEY.EQ.1 .OR. NKEY.EQ.2) THEN
C
C        ELEMENT OR SPECIES DATA
C
         ITASK = NKEY
         IF (NSUB .EQ. 1) GO TO 100
C
         DO 25 N = 2, NSUB
            SUB(N-1) = ' '
            SUB(N-1) = SUB(N)
   25    CONTINUE
         NSUB = NSUB-1
C
      ELSEIF (NKEY .EQ. 3) THEN
C
C        THERMODYNAMIC DATA
C
         IF (NSUB .GT. 1) THEN
            IF ( UPCASE(SUB(2), 3) .EQ. 'ALL') THEN
               THERMO = .FALSE.
               READ (LIN,'(A)') LINE
               CALL IPPARR (LINE, -1, 3, VALUE, NVAL, IER, LOUT)
               IF (NVAL .NE. 3 .OR. IER.NE.0) THEN
                  KERR = .TRUE.
                  WRITE (LOUT, 333)
               ELSE
                  TLO = VALUE(1)
                  TMID = VALUE(2)
                  THI = VALUE(3)
               ENDIF
            ENDIF
         ELSE
C
C           USE THERMODYNAMIC DATABASE FOR DEFAULT TLO,TMID,THI
            OPEN (LTHRM, FORM='FORMATTED', STATUS='UNKNOWN', 
     1                   FILE='therm.dat')
C
            READ (LTHRM,'(A)') LINE
            READ (LTHRM,'(A)') LINE
            CALL IPPARR (LINE, -1, 3, VALUE, NVAL, IER, LOUT)
            IF (NVAL .NE. 3 .OR. IER.NE.0) THEN
               KERR = .TRUE.
               WRITE (LOUT, 333)
            ELSE
               TLO = VALUE(1)
               TMID = VALUE(2)
               THI = VALUE(3)
            ENDIF
            CLOSE (LTHRM)
         ENDIF
C
         CALL CKTHRM (LIN, MDIM, ENAME, MM, AWT, KNAME, KK, KNCF,
     1                KPHSE, KCHRG, WTM, MAXTP, NT, NTR, TLO, TMID,
     2                THI, T, NPCP2, A, ITHRM, KERR, LOUT, LINE)
C
         IF (.NOT. THERMO)
     1      CALL CKPRNT (MDIM, MAXTP, MM, ENAME, KK, KNAME, WTM, KPHSE,
     2                   KCHRG, NT, T, TLO, TMID, THI, KNCF, ITHRM,
     3                   LOUT, KERR)
         I1 = IFIRCH(LINE)
         IF (UPCASE(LINE(I1:), 4) .EQ. 'REAC') GO TO 105
C
      ELSEIF (NKEY .EQ. 4) THEN
C
         ITASK = 4
C        START OF REACTIONS; ARE UNITS SPECIFIED?
         CALL CKUNIT (LINE(:ILEN), AUNITS, EUNITS, IUNITS)
C
         IF (THERMO) THEN
C
C           THERMODYNAMIC DATA
            OPEN (LTHRM, FORM='FORMATTED', STATUS='UNKNOWN', 
     1                   FILE='therm.dat')
            READ (LTHRM,'(A)') LINE
            READ (LTHRM,'(A)') LINE
            CALL IPPARR (LINE, -1, 3, VALUE, NVAL, IER, LOUT)
            IF (NVAL .NE. 3 .OR. IER.NE.0) THEN
               KERR = .TRUE.
               WRITE (LOUT, 333)
            ELSE
               TLO = VALUE(1)
               TMID = VALUE(2)
               THI = VALUE(3)
            ENDIF
            CALL CKTHRM (LTHRM, MDIM, ENAME, MM, AWT, KNAME, KK, KNCF,
     1                   KPHSE, KCHRG, WTM, MAXTP, NT, NTR, TLO, TMID,
     2                   THI, T, NPCP2, A, ITHRM, KERR, LOUT, LINE)
            CALL CKPRNT (MDIM, MAXTP, MM, ENAME, KK, KNAME, WTM, KPHSE,
     1                   KCHRG, NT, T, TLO, TMID, THI, KNCF, ITHRM,
     2                   LOUT, KERR)
            THERMO = .FALSE.
            CLOSE (LTHRM)
         ENDIF
C
         WRITE (LOUT, 1800)
         GO TO 100
      ENDIF
C
      IF (ITASK .EQ. 1) THEN
C
C        ELEMENT DATA
C
         IF (MM .EQ. 0) THEN
            WRITE (LOUT, 200)
            WRITE (LOUT, 300)
            WRITE (LOUT, 200)
         ENDIF
C
         IF (NSUB .GT. 0) THEN
            M1 = MM +1
            CALL CKCHAR (SUB, NSUB, MDIM, ENAME, AWT, MM, KERR, LOUT)
            DO 110 M = M1, MM
               IF (AWT(M) .LE. 0) CALL CKAWTM (ENAME(M), AWT(M))
               WRITE (LOUT, 400) M,ENAME(M)(:4),AWT(M)
               IF (AWT(M) .LE. 0) THEN
                  KERR = .TRUE.
                  WRITE (LOUT, 1000) ENAME(M)
               ENDIF
  110       CONTINUE
         ENDIF
C
      ELSEIF (ITASK .EQ. 2) THEN
C
C        PROCESS SPECIES DATA
C
         IF (KK .EQ. 0) WRITE (LOUT, 200)
         IF (NSUB .GT. 0)
     1   CALL CKCHAR (SUB, NSUB, KDIM, KNAME, WTM, KK, KERR, LOUT)
C
      ELSEIF (ITASK .EQ. 4) THEN
C
C        PROCESS REACTION DATA
C
         IND = 0
         DO 120 N = 1, NSUB
            IND = MAX(IND, INDEX(SUB(N),'/'))
            IF (UPCASE(SUB(N), 3) .EQ. 'DUP') IND = MAX(IND,1)
  120    CONTINUE
         IF (IND .GT. 0) THEN
C
C           AUXILIARY REACTION DATA
C
            CALL CKAUXL (SUB, NSUB, II, KK, KNAME, LOUT, MAXSP, NPAR,
     1                   NSPEC, NTHB, ITHB, NTBS, MAXTB, NKTB, AIK,
     2                   NFAL, IFAL, IDUP, NFAR, PFAL, IFOP, NLAN,
     3                   ILAN, NLAR, PLAN, NREV, IREV, RPAR, NRLT, IRLT, 
     4                   RLAN, NWL, IWL, WL, KERR, NORD, IORD, MAXORD, 
     5                   KORD, RORD, NUNK, NU, NRNU, IRNU, RNU)
C
         ELSE
C
C           THIS IS A REACTION STRING
C
            IF (II .LT. IDIM) THEN
C
               IF (II .GT. 0)
C
C              CHECK PREVIOUS REACTION FOR COMPLETENESS
C
     1         CALL CPREAC (II, MAXSP, NSPEC, NPAR, PAR, RPAR,
     2                      AUNITS, EUNITS, NREAC, NUNK, NU, KCHRG,
     3                      MDIM, MM, KNCF, IDUP, NFAL, IFAL, KFAL,
     4                      NFAR, PFAL, IFOP, NREV, IREV, NTHB, ITHB,
     5                      NLAN, ILAN, NRLT, IRLT, KERR, LOUT, NRNU,
     6                      IRNU, RNU, CKMIN)
C
C              NEW REACTION
C
               II = II+1
               CALL CKREAC (LINE(:ILEN), II, KK, KNAME, LOUT, MAXSP,
     1                      NSPEC, NREAC, NUNK, NU, NPAR, PAR,
     2                      NTHB, ITHB, NFAL, IFAL, KFAL, NWL,
     3                      IWL, WL, NRNU, IRNU, RNU, KERR)
C
            ELSE
               WRITE (LOUT, 1070)
               KERR = .TRUE.
            ENDIF
C
         ENDIF
      ENDIF
      GO TO 100
C
 5000 CONTINUE
C
C     END OF INPUT
C
      IF (II .GT. 0) THEN
C
C              CHECK FINAL REACTION FOR COMPLETENESS
C
          CALL CPREAC (II, MAXSP, NSPEC, NPAR, PAR, RPAR, AUNITS,
     1                 EUNITS, NREAC, NUNK, NU, KCHRG, MDIM, MM,
     2                 KNCF, IDUP, NFAL, IFAL, KFAL, NFAR, PFAL, IFOP,
     3                 NREV, IREV, NTHB, ITHB, NLAN, ILAN, NRLT,
     4                 IRLT, KERR, LOUT, NRNU, IRNU, RNU, CKMIN)
C
C              CHECK REACTIONS DECLARED AS DUPLICATES
C
         DO 500 I = 1, II
            IF (IDUP(I) .LT. 0) THEN
               KERR = .TRUE.
               WRITE (LOUT, 1095) I
            ENDIF
  500    CONTINUE
C
         WRITE (LOUT, '(/1X,A)') ' NOTE: '//IUNITS(:ILASCH(IUNITS))
C
      ELSEIF (THERMO) THEN
C
C        THERE WAS NO REACTION DATA, MAKE SURE SPECIES DATA IS COMPLETE
         OPEN (LTHRM, FORM='FORMATTED', STATUS='UNKNOWN', 
     1                FILE='therm.dat')
C
         READ (LTHRM,'(A)') LINE
         READ (LTHRM,'(A)') LINE
         CALL IPPARR (LINE, -1, 3, VALUE, NVAL, IER, LOUT)
         IF (NVAL .NE. 3 .OR. IER.NE.0) THEN
            KERR = .TRUE.
            WRITE (LOUT, 333)
         ELSE
            TLO = VALUE(1)
            TMID = VALUE(2)
            THI = VALUE(3)
         ENDIF
         CALL CKTHRM (LTHRM, MDIM, ENAME, MM, AWT, KNAME, KK, KNCF,
     1                KPHSE, KCHRG, WTM, MAXTP, NT, NTR, TLO, TMID,
     2                THI, T, NPCP2, A, ITHRM, KERR, LOUT, LINE)
         CALL CKPRNT (MDIM, MAXTP, MM, ENAME, KK, KNAME, WTM, KPHSE,
     1                KCHRG, NT, T, TLO, TMID, THI, KNCF, ITHRM,
     2                LOUT, KERR)
         CLOSE  (LTHRM)
      ENDIF
C
      IF (KERR) THEN
C
         WRITE (LOUT, '(//A)')
     1   ' WARNING...THERE IS AN ERROR IN THE LINKING FILE'
          DO 1150 K = 1, KK
            IF (KCHRG(K) .NE. 0) NCHRG = NCHRG+1
 1150    CONTINUE
         STOP
      ENDIF
C
      LENICK = 1 + (3 + MM)*KK + (2 + 2*MAXSP)*II + NLAN + NRLT
     1           + 3*NFAL + (2 + MAXTB)*NTHB + NREV + NWL + NRNU
     2           + NORD*(1 + MAXORD)
      LENCCK = MM + KK
      LENRCK = 3 + MM + KK*(5 + MAXTP + NTR*NPCP2) + II*7 + NREV
     1           + NPAR*(II + NREV) + NLAR*(NLAN + NRLT)
     2           + NFAR*NFAL + MAXTB*NTHB + NWL + NRNU*MAXSP
     3           + NORD*MAXORD
C
C     OPEN LINKING FILE
C
      OPEN (LINC, FORM='UNFORMATTED', STATUS='UNKNOWN', 
     1            FILE='chem.bin')
C
      WRITE (LINC) VERS, PREC, KERR
      WRITE (LINC) LENICK, LENRCK, LENCCK, MM, KK, II, MAXSP,
     1             MAXTB, MAXTP, NPC, NPAR, NLAR, NFAR, NREV, NFAL,
     2             NTHB, NLAN, NRLT, NWL, NCHRG, NRNU, NORD,
     3             MAXORD, CKMIN
      WRITE (LINC) (ENAME(M), AWT(M), M = 1, MM)
      WRITE (LINC) (KNAME(K), (KNCF(M,K),M=1,MM), KPHSE(K),
     1              KCHRG(K), WTM(K), NT(K), (T(L,K),L=1,MAXTP),
     2              ((A(M,L,K), M=1,NPCP2), L=1,NTR), K = 1, KK)
C
      IF (II .GT. 0) THEN
C
         WRITE (LINC) (NSPEC(I), NREAC(I), (PAR(N,I), N = 1, NPAR),
     1         (NU(M,I), NUNK(M,I), M = 1, MAXSP), I = 1, II)
C
         IF (NREV .GT. 0) WRITE (LINC)
     1      (IREV(N),(RPAR(L,N),L=1,NPAR),N=1,NREV)
C
         IF (NFAL .GT. 0) WRITE (LINC)
     1      (IFAL(N),IFOP(N),KFAL(N),(PFAL(L,N),L=1,NFAR), N = 1, NFAL)
C
         IF (NTHB .GT. 0) WRITE (LINC)
     1      (ITHB(N),NTBS(N),(NKTB(M,N),AIK(M,N),M=1,MAXTB),N=1,NTHB)
C
         IF (NLAN .GT. 0) WRITE (LINC)
     1      (ILAN(N), (PLAN(L,N), L = 1, NLAR), N = 1, NLAN)
C
         IF (NRLT .GT. 0) WRITE (LINC)
     1      (IRLT(N), (RLAN(L,N), L = 1, NLAR), N=1,NRLT)
C
         IF (NWL .GT. 0) WRITE (LINC) (IWL(N), WL(N), N = 1, NWL)
C
         IF (NRNU .GT. 0) WRITE (LINC)
C
C            NRNU, total number of reactions with real stochio. coeff.
C
     1      (IRNU(N), (RNU(M,N), M = 1, MAXSP), N = 1, NRNU)
C
C            IRNU, indices of reaction numbers
C            RNU,   matrix of real stochiometric coefficients
C
         IF (NORD .GT. 0) WRITE (LINC)
C
C            NORD, total number of reactions which use "ORDER"
C
     1      (IORD(N), (KORD(L,N), RORD(L,N), L=1, MAXORD), N=1,NORD)
C
C            IORD, array of reaction numbers
C            KORD, array of species numbers with "ORDER" specified,
C                  -K for forward species, K for reverse species 
C            RORD, array of order coefficients
      ELSE
         WRITE (LOUT, '(/A)')
     1      ' WARNING...NO REACTION INPUT FOUND; ',
     2      ' LINKING FILE HAS NO REACTION INFORMATION ON IT.'
      ENDIF
C
      WRITE (LOUT, '(///A)')
     1   ' NO ERRORS FOUND ON INPUT...CHEMKIN LINKING FILE WRITTEN.'
C
      WRITE (LOUT, '(/A,3(/A,I6))')
     1      ' WORKING SPACE REQUIREMENTS ARE',
     2      '    INTEGER:   ',LENICK,
     3      '    REAL:      ',LENRCK,
     4      '    CHARACTER: ',LENCCK
      CLOSE (LINC)
      CLOSE (LIN)
      CLOSE (LOUT)
C
C----------------------------------------------------------------------C
C
C     FORMATS
C
  200 FORMAT (26X,20('-'))
  300 FORMAT (26X,'ELEMENTS',5X,'ATOMIC',/26X,'CONSIDERED',3X,'WEIGHT')
  333 FORMAT (/6X,'Error...no TLO,TMID,THI given for THERMO ALL...'/)
  400 FORMAT (25X,I3,'. ',A4,G15.6)
C
 1000 FORMAT (6X,'Error...no atomic weight for element ',A)
 1070 FORMAT (6X,'Error...more than IDIM reactions...')
 1095 FORMAT (6X,'Error...no duplicate declared for reaction no.',I3)
 1800 FORMAT (///54X, '(k = A T**b exp(-E/RT))',/,
     1        6X,'REACTIONS CONSIDERED',30X,'A',8X,'b',8X,'E',/)
C
      STOP
      END
C----------------------------------------------------------------------C
      SUBROUTINE CKCHAR (SUB, NSUB, NDIM, STRAY, RAY, NN, KERR, LOUT)
C
C     Extracts names and real values from an array of CHAR*(*)
C     substrings; stores names in STRAY array, real values in RAY;
C     i.e. can be used to store element and atomic weight data,
C     species names, etc.
C
C     Input:   SUB(N),N=1,NSUB  - array of CHAR*(*) substrings
C              NSUB             - number of substrings
C              NDIM             - size of STRAY,RAY arrays
C              NN               - actual number of STRAY found
C              STRAY(N),N=1,NN  - CHAR*(*) array
C              RAY(N),N=1,NN    - Real array
C              LOUT             - output unit for error messages
C     Output:  NN               - incremented if more STRAY found
C              STRAY(N),N=1,NN  - incremented array of STRAY
C              RAY(N),N=1,NN    - incremented array of reals
C              KERR             - logical, .TRUE. = error in data
C
C                                       F. Rupley, Div. 8245, 2/5/88
C----------------------------------------------------------------------C
C*****precision > double
       IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION RAY(*), PAR(1)
      CHARACTER SUB(*)*(*), STRAY(*)*(*), ISTR*80, UPCASE*4
      LOGICAL KERR
C
      ILEN = LEN(STRAY(1))
C
      DO 200 N = 1, NSUB
         IF ( UPCASE(SUB(N), 3) .EQ. 'END') RETURN
         ISTR = ' '
         I1 = INDEX(SUB(N),'/')
         IF (I1 .EQ .1) THEN
            KERR = .TRUE.
            WRITE (LOUT, 130) SUB(N)(:ILASCH(SUB(N)))
         ELSE
            IF (I1 .LE. 0) THEN
               ISTR = SUB(N)
            ELSE
               ISTR = SUB(N)(:I1-1)
            ENDIF
            CALL CKCOMP (ISTR, STRAY, NN, INUM)
C
            IF (INUM .GT. 0) THEN
               WRITE (LOUT, 100) SUB(N)(:ILASCH(SUB(N)))
            ELSE
               IF (NN .LT. NDIM) THEN
                  IF (ISTR(ILEN+1:) .NE. ' ') THEN
                     WRITE (LOUT, 120) SUB(N)(:ILASCH(SUB(N)))
                     KERR = .TRUE.
                  ELSE
                     NN = NN + 1
                     STRAY(NN) = ' '
                     STRAY(NN) = ISTR(:ILEN)
                     IF (I1 .GT. 0) THEN
                        I2 = I1 + INDEX(SUB(N)(I1+1:),'/')
                        ISTR = ' '
                        ISTR = SUB(N)(I1+1:I2-1)
                        CALL IPPARR (ISTR, 1, 1, PAR, NVAL, IER, LOUT)
                        IF (IER .EQ. 0) THEN
                           RAY(NN) = PAR(1)
                        ELSE
                           KERR = .TRUE.
                        ENDIF
                     ENDIF
                  ENDIF
               ELSE
                  WRITE (LOUT, 110) SUB(N)(:ILASCH(SUB(N)))
                  KERR = .TRUE.
               ENDIF
            ENDIF
         ENDIF
  200 CONTINUE
C
  100 FORMAT (6X,'Warning...duplicate array element ignored...',A)
  110 FORMAT (6X,'Error...character array size too small for  ...',A)
  120 FORMAT (6X,'Error...character array element name too long...',A)
  130 FORMAT (6X,'Error...misplaced value...',A)
      END
C----------------------------------------------------------------------C
      SUBROUTINE CKAWTM (ENAME, AWT)
C
C     Returns atomic weight of element ENAME.
C     Input:   ENAME - CHAR*(*) element name
C     Output:  AWT   - real atomic weight
C
C                                       F. Rupley, Div. 8245, 11/11/86
C----------------------------------------------------------------------C
C*****precision > double
       IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      PARAMETER (NATOM = 102)
      DIMENSION ATOM(NATOM)
      CHARACTER ENAME*(*), IATOM(NATOM)*2, UPCASE*2
C
      DATA (IATOM(I),ATOM(I),I=1,40) /
     *'H ',  1.00797, 'HE',  4.00260, 'LI',  6.93900, 'BE',  9.01220,
     *'B ', 10.81100, 'C ', 12.01115, 'N ', 14.00670, 'O ', 15.99940,
     *'F ', 18.99840, 'NE', 20.18300, 'NA', 22.98980, 'MG', 24.31200,
     *'AL', 26.98150, 'SI', 28.08600, 'P ', 30.97380, 'S ', 32.06400,
     *'CL', 35.45300, 'AR', 39.94800, 'K ', 39.10200, 'CA', 40.08000,
     *'SC', 44.95600, 'TI', 47.90000, 'V ', 50.94200, 'CR', 51.99600,
     *'MN', 54.93800, 'FE', 55.84700, 'CO', 58.93320, 'NI', 58.71000,
     *'CU', 63.54000, 'ZN', 65.37000, 'GA', 69.72000, 'GE', 72.59000,
     *'AS', 74.92160, 'SE', 78.96000, 'BR', 79.90090, 'KR', 83.80000,
     *'RB', 85.47000, 'SR', 87.62000, 'Y ', 88.90500, 'ZR', 91.22000/
C
      DATA (IATOM(I),ATOM(I),I=41,80) /
     *'NB', 92.90600, 'MO', 95.94000, 'TC', 99.00000, 'RU',101.07000,
     *'RH',102.90500, 'PD',106.40000, 'AG',107.87000, 'CD',112.40000,
     *'IN',114.82000, 'SN',118.69000, 'SB',121.75000, 'TE',127.60000,
     *'I ',126.90440, 'XE',131.30000, 'CS',132.90500, 'BA',137.34000,
     *'LA',138.91000, 'CE',140.12000, 'PR',140.90700, 'ND',144.24000,
     *'PM',145.00000, 'SM',150.35000, 'EU',151.96000, 'GD',157.25000,
     *'TB',158.92400, 'DY',162.50000, 'HO',164.93000, 'ER',167.26000,
     *'TM',168.93400, 'YB',173.04000, 'LU',174.99700, 'HF',178.49000,
     *'TA',180.94800, 'W ',183.85000, 'RE',186.20000, 'OS',190.20000,
     *'IR',192.20000, 'PT',195.09000, 'AU',196.96700, 'HG',200.59000/
C
      DATA (IATOM(I),ATOM(I),I=81,NATOM) /
     *'TL',204.37000, 'PB',207.19000, 'BI',208.98000, 'PO',210.00000,
     *'AT',210.00000, 'RN',222.00000, 'FR',223.00000, 'RA',226.00000,
     *'AC',227.00000, 'TH',232.03800, 'PA',231.00000, 'U ',238.03000,
     *'NP',237.00000, 'PU',242.00000, 'AM',243.00000, 'CM',247.00000,
     *'BK',249.00000, 'CF',251.00000, 'ES',254.00000, 'FM',253.00000,
     *'D ',002.01410, 'E',5.45E-4/
C
      CALL CKCOMP ( UPCASE(ENAME, 2), IATOM, NATOM, L)
      IF (L .GT. 0) AWT = ATOM(L)
      RETURN
      END
C----------------------------------------------------------------------C
      SUBROUTINE CKTHRM (LUNIT, MDIM, ENAME, MM, AWT, KNAME, KK, KNCF,
     1                   KPHSE, KCHRG, WTM, MAXTP, NT, NTR, TLO, TMID,
     2                   THI, T, NPCP2, A, ITHRM, KERR, LOUT, ISTR)
C
C     Finds thermodynamic data and elemental composition for species
C     Input:  LUNIT  - unit number for input of thermo properties
C             MDIM   - maximum number of elements allowed
C             ENAME(M),M=1,MM  - array of CHAR*(*) element names
C             MM     - total number of elements declared
C             AWT(M),M=1,MM    - array of atomic weights for elements
C             KNAME(K),K=1,KK  - array of CHAR*(*) species names
C             KK     - total number of species declared
C             LOUT   - output unit for messages
C             NT(K),K=1,KK - number of temperature values
C             NTR - number of temperature ranges
C     Output: KNCF(M,K) - elemental composition of species
C             KPHSE(K),K=1,KK - integer array, species phase
C             KCHRG(K),K=1,KK - integer array of species charge
C                      =0, if no electrons,
C                      =(-1)*number of electrons present
C             WTM(K),K=1,KK - array of molecular weights of species
C             A(M,L,K)- array of thermodynamic coefficients
C             T(N),N=1,NT - array of temperatures
C             KERR   - logical error flag
C----------------------------------------------------------------------C
C*****precision > double
       IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION WTM(*), NT(*), T(MAXTP,*), KPHSE(*), KNCF(MDIM,*),
     1          KCHRG(*), A(NPCP2,NTR,*), AWT(*), VALUE(5)
      CHARACTER ENAME(*)*(*), KNAME(*)*(*), LINE(4)*80, ELEM*16
      CHARACTER UPCASE*4, ISTR*80, SUB(80)*80
      LOGICAL KERR, ITHRM(*)
C
      IF (MM.LE.0 .OR. KK.LE.0) WRITE (LOUT, 80)
C
      GO TO 20
   10 CONTINUE
      ISTR = ' '
      READ (LUNIT,'(A)',END=40) ISTR
   20 CONTINUE
      ILEN = IPPLEN(ISTR)
      IF (ILEN .LE. 0) GO TO 10
C
      CALL CKISUB (ISTR(:ILEN), SUB, NSUB)
      CALL CKCOMP (SUB(1), KNAME, KK, K)
      IF (K .EQ. 0) THEN
         IF (UPCASE(SUB(1), 3) .EQ. 'END' .OR.
     1       UPCASE(SUB(1), 4) .EQ. 'REAC') RETURN
         GO TO 10
      ENDIF
C
      IF (ITHRM(K)) GO TO 10
      ITHRM(K) = .TRUE.
      LINE(1) = ' '
      LINE(1) = ISTR
      DO 25 L = 2, 4
         LINE(L) = ' '
         READ (LUNIT,'(A)',END=40) LINE(L)
   25 CONTINUE
C
      ICOL = 20
      DO 60 I = 1, 5
         ICOL = ICOL + 5
         IF (I .EQ. 5) ICOL = 74
         ELEM  = LINE(1)(ICOL:ICOL+1)
         IELEM = 0
C
         IF (LINE(1)(ICOL+2:ICOL+4) .NE. ' ') THEN
            CALL IPPARR
     1      (LINE(1)(ICOL+2:ICOL+4), 0, 1, VALUE, NVAL, IER, LOUT)
            IELEM = VALUE(1)
         ENDIF
C
         IF (ELEM.NE.' ' .AND. IELEM.NE.0) THEN
            IF (UPCASE(ELEM, 1) .EQ. 'E')
     1             KCHRG(K)=KCHRG(K)+IELEM*(-1)
            CALL CKCOMP (ELEM, ENAME, MM, M)
            IF (M .GT. 0) THEN
               KNCF(M,K) = IELEM
               WTM(K) = WTM(K) + AWT(M)*FLOAT(IELEM)
            ELSE
               WRITE (LOUT, 100) ELEM,KNAME(K)(:10)
               KERR = .TRUE.
            ENDIF
         ENDIF
   60 CONTINUE
C
      IF (UPCASE(LINE(1)(45:),1) .EQ. 'L') KPHSE(K)=1
      IF (UPCASE(LINE(1)(45:),1) .EQ. 'S') KPHSE(K)=-1
C
C-----Currently allows for three temperatures, two ranges;
C     in future, NT(K) may vary, NTR = NT(K)-1
C
      T(1,K) = TLO
      IF (LINE(1)(46:55) .NE. ' ') CALL IPPARR
     1   (LINE(1)(46:55), 0, 1, T(1,K), NVAL, IER, LOUT)
C
      T(2,K) = TMID
      IF (LINE(1)(66:73) .NE. ' ') CALL IPPARR
     1   (LINE(1)(66:73), 0, 1, T(2,K), NVAL, IER, LOUT)
C
      T(NT(K),K) = THI
      IF (LINE(1)(56:65) .NE. ' ') CALL IPPARR
     1   (LINE(1)(56:65), 0, 1, T(NT(K),K), NVAL, IER, LOUT)
C
      READ (LINE(2)(:75),'(5E15.8)') (A(I,NTR,K),I=1,5)
      READ (LINE(3)(:75),'(5E15.8)')
     1            (A(I,NTR,K),I=6,7),(A(I,1,K),I=1,3)
      READ (LINE(4)(:60),'(4E15.8)') (A(I,1,K),I=4,7)
      GO TO 10
C
   40 RETURN
   80 FORMAT (6X,'Warning...THERMO cards misplaced will be ignored...')
  100 FORMAT (6X,'Error...element...',A,'not declared for...',A)
      END
C----------------------------------------------------------------------C
      SUBROUTINE CKREAC (LINE, II, KK, KNAME, LOUT, MAXSP, NSPEC, NREAC,
     1                   NUNK, NU, NPAR, PAR, NTHB, ITHB,
     2                   NFAL, IFAL, KFAL, NWL, IWL, WL, 
     3                   NRNU, IRNU, RNU, KERR)
C
C     CKREAC parses the main CHAR*(*) line representing a gas-phase
C     reaction; first, the real Arrhenius parameters are located and
C     stored in PAR(N,I),N=1,NPAR, where I is the reaction number;
C     then a search is made over the reaction string:
C
C     '=','<=>': reaction I is reversible;
C     '=>'     : reaction I is irreversible;
C
C     '(+[n]KNAME(K))': reaction I is a fall-off reaction;
C                       NFAL is incremented, the total number of
C                       fall-off reactions;
C                       IFAL(NFAL)=I, KFAL(NFAL)=K;
C                       this species is eliminated from consideration
C                       as a reactant or product in this reaction.
C
C     '(+M)'   : reaction I is a fall-off reaction;
C                NFAL is incremented, IFAL(NFAL)=I, KFAL(NFAL)=0;
C
C     '+[n]KNAME(K)': NSPEC(I) is incremented, the total number of
C                     species for this reaction;
C                     n is an optional stoichiometric coefficient
C                     of KNAME(K), if omitted, n=1;
C                     if this string occurs before the =/-,
C                     NREAC(I) is incremented, the total number of
C                     reactants for this reaction, NUNK(N,I)=K, and
C                     NU(N,I) = -n, where N=1-3 is reserved for
C                     reactants;
C                     if this string occurs after the =/-,
C                     NUNK(N,I) = K, and NU(N,I) = n, where N=4-6
C                     is reserved for products;
C
C     '+M' : I is a third-body reaction; NTHB is incremented, the
C            total number of third-body reactions, and ITHB(NTHB)=I.
C
C     Input:  LINE  - a CHAR*(*) line (from data file)
C             II    - the index of this reaction, and the total number
C                     of reactions found so far.
C             KK    - actual integer number of species
C             KNAME(K),K=1,KK - array of CHAR*(*) species names
C             LOUT  - output unit for error messages
C             MAXSP - maximum number of species allowed in reaction
C             NPAR  - number of parameters expected
C     A '!' will comment out a line, or remainder of the line.
C
C     Output: NSPEC - total number of reactants+products in reaction
C             NREAC - number of reactants
C             NUNK  - species numbers for the NSPEC species
C             NU    - stoichiometric coefficients for the NSPEC spec.
C             NFAL  - total number of fall-off reactions
C             IFAL  - reaction numbers for the NFAL reactions
C             KFAL  - 3rd body species numbers for the NFAL reactions
C             NTHB  - total number of 3rd-body reactions
C             ITHB  - reaction numbers for the NTHB reactions
C             NWL   - number of radiation-enhanced reactions
C             IWL   - reaction numbers for the NWL reactions
C             WL    - radiation wavelengths for the NWL reactions
C             KERR  - logical, .TRUE. = error in data file
C
C                                      F. Rupley, Div. 8245, 5/13/86
C----------------------------------------------------------------------C
C*****precision > double
       IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION NSPEC(*), NREAC(*), NUNK(MAXSP,*), NU(MAXSP,*),
     1          PAR(NPAR,*), IFAL(*), KFAL(*), ITHB(*), IWL(*), WL(*),
     2          IRNU(*), RNU(MAXSP,*), IPLUS(20)
      CHARACTER KNAME(*)*(*), LINE*(*), CNUM(11)*1, UPCASE*4
      CHARACTER*80 ISTR, IREAC, IPROD, ISPEC, INAME, ITEMP
      LOGICAL KERR, LTHB, LWL, LRSTO
      DATA CNUM/'.','0','1','2','3','4','5','6','7','8','9'/
C
      LTHB = .FALSE.
      LWL = .FALSE.
      NSPEC(II) = 0
      NREAC(II) = 0
C
C----------Find NPAR real parameters------------------------
C
      CALL IPNPAR (LINE, NPAR, ISTR, ISTART)
      CALL IPPARR (ISTR, 1, NPAR, PAR(1,II), NVAL, IER, LOUT)
      IF (IER .NE. 0) KERR = .TRUE.
C
C-----Remove blanks from reaction string
C
      INAME = ' '
      ILEN = 0
      DO 10 I = 1, ISTART-1
         IF (LINE(I:I) .NE. ' ') THEN
            ILEN = ILEN+1
            INAME(ILEN:ILEN) = LINE(I:I)
         ENDIF
   10 CONTINUE
C
C-----Find reaction string, product string
C
      I1 = 0
      I2 = 0
      DO 25 I = 1, ILEN
         IF (I1 .LE. 0) THEN
            IF (INAME(I:I+2) .EQ. '<=>') THEN
               I1 = I
               I2 = I+2
               IR = 1
            ELSEIF (INAME(I:I+1) .EQ. '=>') THEN
               I1 = I
               I2 = I+1
               IR = -1
            ELSEIF (I.GT.1 .AND. INAME(I:I).EQ.'='
     1                  .AND. INAME(I-1:I-1).NE.'=') THEN
               I1 = I
               I2 = I
               IR = 1
            ENDIF
         ENDIF
   25 CONTINUE
C
      IF (ILASCH(INAME).GE.45 .AND. I1.GT.0) THEN
         WRITE (LOUT, 1900) II,INAME(:I1-1),(PAR(N,II),N=1,NPAR)
         WRITE (LOUT, 1920) INAME(I1:)
      ELSE
          WRITE (LOUT, 1900) II,INAME(:45),(PAR(N,II),N=1,NPAR)
      ENDIF
C
      IREAC = ' '
      IPROD = ' '
      IF (I1 .GT. 0) THEN
         IREAC = INAME(:I1-1)
         IPROD = INAME(I2+1:)
      ELSE
C
C-----did not find delimiter
C
         WRITE (LOUT, 660)
         KERR = .TRUE.
         RETURN
      ENDIF
C
      LRSTO = ((INDEX(IREAC,'.').GT.0) .OR. (INDEX(IPROD,'.').GT.0))
      IF (LRSTO) THEN
         NRNU = NRNU + 1
         IRNU(NRNU) = II
      ENDIF
C
      IF (INDEX(IREAC,'=>').GT.0 .OR. INDEX(IPROD,'=>').GT.0) THEN
C
C-----more than one '=>'
C
         WRITE (LOUT, 800)
         KERR = .TRUE.
         RETURN
      ENDIF
C
C-----Is this a fall-off reaction?
C
      IF (INDEX(IREAC,'(+').GT.0 .OR. INDEX(IPROD,'(+').GT.0) THEN
         KRTB = 0
         KPTB = 0
         DO 300 J = 1, 2
            ISTR = ' '
            KTB  = 0
            IF (J .EQ. 1) THEN
               ISTR = IREAC
            ELSE
               ISTR = IPROD
            ENDIF
C
            DO 35 N = 1, ILASCH(ISTR)-1
               IF (ISTR(N:N+1) .EQ. '(+') THEN
                  I1 = N+2
                  I2 = I1 + INDEX(ISTR(I1:),')')-1
                  IF (I2 .GT. I1) THEN
                     IF (ISTR(I1:I2-1).EQ.'M' .OR.
     1                   ISTR(I1:I2-1).EQ.'m') THEN
                         IF (KTB .NE. 0) THEN
                            WRITE (LOUT, 630)
                            KERR = .TRUE.
                            RETURN
                         ELSE
                            KTB = -1
                         ENDIF
                     ELSE
                        CALL CKCOMP (ISTR(I1:I2-1), KNAME, KK, KNUM)
                        IF (KNUM .GT. 0) THEN
                           IF (KTB .NE. 0) THEN
                              WRITE (LOUT, 630)
                              KERR = .TRUE.
                              RETURN
                           ELSE
                              KTB = KNUM
                           ENDIF
                        ENDIF
                     ENDIF
                     IF (KTB .NE. 0) THEN
                        ITEMP = ' '
                        IF (I1 .EQ. 1) THEN
                           ITEMP = ISTR(I2+1:)
                        ELSE
                           ITEMP = ISTR(:I1-3)//ISTR(I2+1:)
                        ENDIF
                        IF (J .EQ. 1) THEN
                           IREAC = ' '
                           IREAC = ITEMP
                           KRTB = KTB
                        ELSE
                           IPROD = ' '
                           IPROD = ITEMP
                           KPTB = KTB
                        ENDIF
                     ENDIF
                  ENDIF
               ENDIF
   35       CONTINUE
  300    CONTINUE
C
         IF (KRTB.NE.0 .OR. KPTB.NE.0) THEN
C
C           does product third-body match reactant third-body
C
            IF (KRTB.LE.0 .AND. KPTB.LE.0) THEN
C
               NFAL = NFAL + 1
               IFAL(NFAL) = II
               KFAL(NFAL) = 0
C
               LTHB = .TRUE.
               NTHB = NTHB + 1
               ITHB(NTHB) = II
C
            ELSEIF (KRTB .EQ. KPTB) THEN
               NFAL = NFAL + 1
               IFAL(NFAL) = II
               KFAL(NFAL) = KRTB
C
            ELSE
C
               WRITE (LOUT, 640)
               KERR = .TRUE.
               RETURN
            ENDIF
         ENDIF
      ENDIF
C
C----------Find reactants, products-------------------------
C
      DO 600 J = 1, 2
         ISTR = ' '
         LTHB = .FALSE.
         IF (J .EQ. 1) THEN
            ISTR = IREAC
            NS = 0
         ELSE
            ISTR = IPROD
            NS = 3
         ENDIF
C
C-----------store pointers to '+'-signs
C
         NPLUS = 1
         IPLUS(NPLUS) = 0
         DO 500 L = 2, ILASCH(ISTR)-1
            IF (ISTR(L:L).EQ.'+') THEN
               NPLUS = NPLUS + 1
               IPLUS(NPLUS) = L
            ENDIF
  500    CONTINUE
         NPLUS = NPLUS + 1
         IPLUS(NPLUS) = ILASCH(ISTR)+1
C
         NSTART = 1
  505    CONTINUE
         N1 = NSTART
         DO 510 N = NPLUS, N1, -1
            ISPEC = ' '
            ISPEC = ISTR(IPLUS(N1)+1 : IPLUS(N)-1)
C
            IF (UPCASE(ISPEC, 1).EQ.'M' .AND.
     1               (ISPEC(2:2).EQ.' ' .OR. ISPEC(2:2).EQ.'+')) THEN
               IF (LTHB) THEN
                  WRITE (LOUT, 900)
                  KERR = .TRUE.
                  RETURN
               ELSEIF (NFAL.GT.0 .AND. IFAL(NFAL).EQ.II) THEN
                  WRITE (LOUT, 640)
                  KERR = .TRUE.
                  RETURN
               ELSE
                  LTHB = .TRUE.
                  IF (NTHB.EQ.0 .OR.
     1               (NTHB.GT.0.AND.ITHB(NTHB).NE.II)) THEN
                      NTHB = NTHB + 1
                      ITHB(NTHB) = II
                  ENDIF
                  IF (N .EQ. NPLUS) GO TO 600
                  NSTART = N
                  GO TO 505
               ENDIF
C
            ELSEIF (UPCASE(ISPEC, 2) .EQ. 'HV') THEN
               IF (LWL) THEN
                  WRITE (LOUT, 670)
                  KERR = .TRUE.
                  RETURN
               ELSE
                  LWL = .TRUE.
                  NWL = NWL + 1
                  IWL(NWL) = II
                  WL(NWL) = 1.0
                  IF (J .EQ. 1) WL(NWL) = -1.0
                  IF (N .EQ. NPLUS) GO TO 600
                  NSTART = N
                  GO TO 505
               ENDIF
            ENDIF
C
C-----------does this string start with a number?
C
            IND = 0
            DO 334 L = 1, LEN(ISPEC)
               NTEST = 0
               DO 333 M = 1, 11
                  IF (ISPEC(L:L) .EQ. CNUM(M)) THEN
                     NTEST=M
                     IND = L
                  ENDIF
  333          CONTINUE
               IF (NTEST .EQ. 0) GO TO 335
  334       CONTINUE
  335       CONTINUE
C
            RVAL = 1.0
            IVAL = 1
            IF (IND .GT. 0) THEN
               IF (LRSTO) THEN
                  CALL IPPARR (ISPEC(:IND), 1, 1, RVAL, NVAL, 
     1                         IER, LOUT)
               ELSE
                  CALL IPPARI (ISPEC(:IND), 1, 1, IVAL, NVAL,
     1                        IER, LOUT)
               ENDIF
               IF (IER .EQ. 0) THEN
                  ITEMP = ' '
                  ITEMP = ISPEC(IND+1:)
                  ISPEC = ' '
                  ISPEC = ITEMP
               ELSE
                  KERR = .TRUE.
                  RETURN
               ENDIF
            ENDIF
C
            CALL CKCOMP (ISPEC, KNAME, KK, KNUM)
            IF (KNUM .EQ. 0) THEN
               IF ((N-N1) .GT. 1) GO TO 510
               WRITE (LOUT, 680) ISPEC(:ILASCH(ISPEC))
               KERR = .TRUE.
            ELSE
C
C--------------a species has been found
C
               IF (J .EQ. 1) THEN
                  IVAL = -IVAL
                  RVAL = -RVAL
               ENDIF
C
C--------------increment species coefficient count
C
               NNUM = 0
               IF (LRSTO) THEN
                  DO 110 K = 1, NS
                     IF (KNUM.EQ.NUNK(K,II) .AND.
     1                   RNU(K,NRNU)/RVAL.GT.0) THEN
                         NNUM = K
                         RNU(NNUM,NRNU) = RNU(NNUM,NRNU) + RVAL
                     ENDIF
  110             CONTINUE
               ELSE
                  DO 111 K = 1, NS
                     IF (KNUM.EQ.NUNK(K,II) .AND.
     1                   NU(K,II)/IVAL.GT.0) THEN
                        NNUM=K
                        NU(NNUM,II) = NU(NNUM,II) + IVAL
                     ENDIF
  111             CONTINUE
               ENDIF
C
               IF (NNUM .LE. 0) THEN
C
C-----------------are there too many species?
C
                  IF (J.EQ.1 .AND. NS.EQ.3) THEN
                     WRITE (LOUT, 690)
                     KERR = .TRUE.
                     RETURN
                  ELSEIF (J.EQ.2 .AND. NS.EQ.MAXSP) THEN
                     WRITE (LOUT, 700)
                     KERR = .TRUE.
                     RETURN
                  ELSE
C
C--------------------increment species count
C
                     NS = NS + 1
                     NSPEC(II) = NSPEC(II)+1
                     IF (J .EQ. 1) NREAC(II) = NS
                     NUNK(NS,II) = KNUM
                     IF (LRSTO) THEN
                        RNU(NS,NRNU) = RVAL
                     ELSE
                        NU(NS,II)   = IVAL
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
            IF (N .EQ. NPLUS) GO TO 600
            NSTART = N
            GO TO 505
C
  510    CONTINUE
  600 CONTINUE
C
      NSPEC(II) = IR*NSPEC(II)
C
  630 FORMAT (6X,'Error...more than one fall-off declaration...')
  640 FORMAT (6X,'Error in fall-off declaration...')
  650 FORMAT (6X,'Error...reaction string not found...')
  660 FORMAT (6X,'Error in reaction...')
  670 FORMAT (6X,'Error in HV declaration...')
  680 FORMAT (6X,'Error...undeclared species...',A)
  690 FORMAT (6X,'Error...more than 3 reactants...')
  700 FORMAT (6X,'Error...more than 3 products...')
  800 FORMAT (6X,'Error in reaction delimiter...')
  900 FORMAT (6X,'Error in third-body declaration...')
C 1900 FORMAT (I4,'. ',A,T51,E10.3,F7.3,F11.3)
 1900 FORMAT (I4,'. ', A, T53, 1PE8.2, 2X, 0PF5.1, 2X, F9.1)
 1920 FORMAT (6X,A)
      RETURN
      END
C----------------------------------------------------------------------C
      SUBROUTINE CKAUXL (SUB, NSUB, II, KK, KNAME, LOUT, MAXSP, NPAR,
     1                   NSPEC, NTHB, ITHB, NTBS, MAXTB, NKTB, AIK, 
     2                   NFAL, IFAL, IDUP, NFAR, PFAL, IFOP, NLAN,
     3                   ILAN, NLAR, PLAN, NREV, IREV, RPAR, NRLT, IRLT, 
     4                   RLAN, NWL, IWL, WL, KERR, NORD, IORD, MAXORD, 
     5                   KORD, RORD, NUNK, NU, NRNU, IRNU, RNU)
C
C     CKAUXL parses the auxiliary CHAR*(*) lines representing
C     additional options for a gas-phase reaction; data is stored
C     based on finding a 'keyword' followed by its required
C     parameters:
C
C     KNAME(K)/val1/: this is an enhanced third-body;
C
C        if ITHB(NTHB) <> I, this is an error, reaction I is not a
C                            third-body reaction;
C        else NTBS(NTHB) is incremented,
C             AIK(NTBS(NTHB),NTHB) = K,
C             NKTB(NTBS(NTHB)),NTHB) = val1;
C
C     (LOW,TROE, and SRI define fall-off data):
C
C     LOW/val1 val2 val3/: PFAL(N,NFAL) = val(N),N=1,3;
C
C        if IFAL(NFAL)<>I, this is an error, reaction I is not a
C                          fall-off reaction;
C        if ILAN(NLAN)=I, this is an error, cannot have T-L numbers.
C        if IRLT(NRLT)=I, this is an error,         "
C        if IREV(NREV)=I, this is an error, cannot declare reverse
C                         parameters;
C        if IFOP(NFAL)>0, this is an error, LOW already declared;
C        else
C           IFOP(NFAL) = ABS(IFOP(NFAL))
C
C     TROE/val1 val2 val3 [val4]/: PFAL(N,NFAL) = val(N),N=4,7;
C
C        if IFAL(NFAL)<>I, this is an error, reaction I is not a
C                          fall-off reaction;
C        if ILAN(NLAN)=I, this is an error, cannot have T-L numbers.
C        if IRLT(NRLT)=I, this is an error,         "
C        if IREV(NREV)=I, this is an error, cannot declare reverse
C                         parameters;
C        if ABS(IFOP(NFAL)).GT.1, this is an error,
C        else
C        if 3 TROE values, IFOP(NFAL) = 3*IFOP(NFAL);
C        if 4 TROE values, IFOP(NFAL) = 4*IFOP(NFAL);
C
C     SRI/val1 val2 val3/: PFAL(N,NFAL) = val(N),N=4,6;
C
C        if IFAL(NFAL)<>I, this is an error, reaction I is not a
C                          fall-off reaction;
C        if ILAN(NLAN)=I, this is an error, cannot have T-L numbers.
C        if IRLT(NRLT)=I, this is an error,         "
C        if IREV(NREV)=I, this is an error, cannot declare reverse
C                         parameters;
C        if ABS(IFOP(NFAL))>1, this is an error;
C        else
C        if IFOP(NFAL)= 2*IFOP(NFAL);
C
C     LT/val1 val2/:
C        if IFAL(NFAL)=I, this is an error, cannot have fall-off and
C                         T-L numbers;
C        else increment NLAN, the number of T-L reactions,
C             ILAN(NLAN)=I, PLAN(N,NLAN)=val(N),N=1,2
C        if IREV(NREV)=I, need IRLT(NRLT)=I.
C
C     REV[ERSE]/val1 val2 val3/ :
C        if IFAL(NFAL)=I, this is an error;
C        if IREV(NREV)=I, this is an error, REV already declared;
C        if NSPEC(I)<0, this an error, as I is irreversible;
C        else increment NREV, the number of reactions with reverse
C             parameters given,
C             IREV(NREV)=I, RPAR(N,NREV)=val(N),N=1,3;
C             if ILAN(NLAN)=I, need IRLT(NRLT)=I;
C             if IRLT(NRLT)=I, need ILAN(NRLT)=I.
C
C     RLT/val1 val2/:
C       if IFAL(NFAL)=I, this is an error, cannot have fall-off and
C                        T-L numbers;
C       if IRLT(NRLT)=I, this is an error, RLT already declared;
C       else increment NRLT, the number of reactions with BOTH
C                      reverse parameters given, and T-L numbers;
C            IRLT(NRLT)=I, RLAN(N,NRLT)=val(N),N=1,2;
C            if IREV(NREV)<>I, need IREV(NREV)=I;
C            if ILAN(NREV)<>I, need ILAN(NLAN)=I;
C
C    DUP[LICATE]:
C       This reaction is allowed to be duplicated.
C
C     Input:  LINE - CHAR*(*) auxiliary information string
C             KK   - total number of species declared
C             KNAME- CHAR*(*) species names
C             LOUT - output unit for error messages
C             MAXSP- maximum third bodies allowed in a reaction
C     Output: NTHB - total number of reactions with third bodies
C             ITHB - array of third-body reaction numbers
C             AIK  - non-zero third body enhancement factors
C             NKTB - array of species numbers for the third body
C                         enchancement factors
C             NFAL - total number of fall-off reactions
C             IFAL - array of fall-off reaction numbers
C             IFOP - array of fall-off type
C             PFAL - fall-off parameters
C             NLAN - total number of Landau-Teller reactions
C             ILAN - array of T-L reaction numbers
C             NLAR - number of Landau-Teller numbers allowed
C             PLAN - array of Landau-Teller numbers
C             NRLT - total number of 'reverse' T-L reactions
C             IRLT - array of 'reverse' T-L reaction numbers
C             RLAN - array of 'reverse' Landau-Teller numbers
C             NWL  - total number of radiation-enhanced reactions
C             IWL  - array of radiation-enhanced reaction numbers
C             WL   - array of wavelengths
C             KERR - logical, = .TRUE. if error found
C                                        F. Rupley, Div. 8245, 5/27/87
C----------------------------------------------------------------------C
C*****precision > double
       IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION NSPEC(*), ITHB(*), NTBS(*), NKTB(MAXTB,*), IDUP(*),
     1          AIK(MAXTB,*), IFAL(*), IFOP(*), PFAL(NFAR,*),
     2          ILAN(*), PLAN(NLAR,*), IREV(*), RPAR(NPAR,*), IRLT(*),
     3          RLAN(NLAR,*), IWL(*), WL(*), VAL(1), IORD(*), 
     4          KORD(MAXORD,*), RORD(MAXORD,*), NUNK(MAXSP,*),
     5          NU(MAXSP,*), IRNU(*), RNU(MAXSP,*)
      CHARACTER SUB(*)*(*), KNAME(*)*(*), KEY*80, RSTR*80, UPCASE*4,
     1          ISTR*80
      LOGICAL KERR, LLAN, LRLT, LTHB, LFAL, LTRO, LSRI, LWL, LREV,
     1        LFORD, LRORD
C
      LTHB = (NTHB.GT.0 .AND. ITHB(NTHB).EQ.II)
      LFAL = (NFAL.GT.0 .AND. IFAL(NFAL).EQ.II)
      LWL  = (NWL .GT.0 .AND. IWL(NWL)  .EQ.II)
      LREV = (NREV.GT.0 .AND. IREV(NREV).EQ.II)
      LLAN = (NLAN.GT.0 .AND. ILAN(NLAN).EQ.II)
      LRLT = (NRLT.GT.0 .AND. IRLT(NRLT).EQ.II)
      LTRO = (NFAL.GT.0 .AND. IFAL(NFAL).EQ.II .AND. IFOP(NFAL).GT.2)
      LSRI = (NFAL.GT.0 .AND. IFAL(NFAL).EQ.II .AND. IFOP(NFAL).EQ.2)
C
      DO 500 N = 1, NSUB
         ILEN = ILASCH(SUB(N))
         KEY = ' '
C
         IF ( UPCASE(SUB(N), 3) .EQ. 'DUP') THEN
            IDUP(II) = -1
            WRITE (LOUT, 4000)
            GO TO 500
         ELSE
            I1 = INDEX(SUB(N),'/')
            I2 = INDEX(SUB(N)(I1+1:),'/')
            IF (I1.LE.0 .OR. I2.LE.0) THEN
               KERR = .TRUE.
               WRITE (LOUT, 2090) SUB(N)(:ILEN)
               GO TO 500
            ENDIF
            KEY = SUB(N)(:I1-1)
            RSTR = ' '
            RSTR = SUB(N)(I1+1:I1+I2-1)
         ENDIF
C
         IF (UPCASE(KEY, 3).EQ.'LOW' .OR.
     1       UPCASE(KEY, 4).EQ.'TROE'.OR.
     2       UPCASE(KEY, 3).EQ.'SRI') THEN
C
C        FALL-OFF DATA
C
            IF ((.NOT.LFAL) .OR. LLAN .OR. LRLT .OR. LREV) THEN
               KERR = .TRUE.
               IF (.NOT. LFAL) WRITE (LOUT, 1050) SUB(N)(:ILEN)
               IF (LLAN)       WRITE (LOUT, 1060) SUB(N)(:ILEN)
               IF (LRLT)       WRITE (LOUT, 1070) SUB(N)(:ILEN)
               IF (LREV)       WRITE (LOUT, 1090) SUB(N)(:ILEN)
            ELSE
C
               IF (UPCASE(KEY, 3) .EQ. 'LOW') THEN
                  IF (IFOP(NFAL) .GT. 0) THEN
                     WRITE (LOUT, 2000) SUB(N)(:ILEN)
                     KERR = .TRUE.
                  ELSE
                     IFOP(NFAL) = ABS(IFOP(NFAL))
                     CALL IPPARR (RSTR,1,3,PFAL(1,NFAL),NVAL,IER,LOUT)
                     IF (IER .NE. 0) KERR = .TRUE.
                     WRITE (LOUT, 3050) (PFAL(L,NFAL),L=1,3)
                  ENDIF
C
               ELSEIF (UPCASE(KEY, 4) .EQ. 'TROE') THEN
                  IF (LTRO .OR. LSRI) THEN
                     KERR = .TRUE.
                     IF (LTRO) WRITE (LOUT, 2010) SUB(N)(:ILEN)
                     IF (LSRI) WRITE (LOUT, 2030) SUB(N)(:ILEN)
                  ELSE
                     LTRO = .TRUE.
                     CALL IPPARR (RSTR,1,-4,PFAL(4,NFAL),NVAL,IER,LOUT)
                     IF (NVAL .EQ. 3) THEN
                        IFOP(NFAL) = 3*IFOP(NFAL)
                        WRITE (LOUT, 3080) (PFAL(L,NFAL),L=4,6)
                     ELSEIF (NVAL .EQ. 4) THEN
                        IFOP(NFAL) = 4*IFOP(NFAL)
                        WRITE (LOUT, 3090) (PFAL(L,NFAL),L=4,7)
                     ELSE
                        WRITE (LOUT, 2020) SUB(N)(:ILEN)
                        KERR = .TRUE.
                     ENDIF
                  ENDIF
C
               ELSEIF (UPCASE(KEY, 3) .EQ. 'SRI') THEN
                  IF (LTRO .OR. LSRI) THEN
                     KERR = .TRUE.
                     IF (LTRO) WRITE (LOUT, 2030) SUB(N)(:ILEN)
                     IF (LSRI) WRITE (LOUT, 2040) SUB(N)(:ILEN)
                  ELSE
                     LSRI = .TRUE.
                     IFOP(NFAL) = 2*IFOP(NFAL)
                     CALL IPPARR (RSTR,1,-5,PFAL(4,NFAL),NVAL,IER,LOUT)
                     IF (NVAL .EQ. 3) THEN
                        PFAL(7,NFAL) = 1.0
                        PFAL(8,NFAL) = 0.0
                        WRITE (LOUT, 3060) (PFAL(L,NFAL),L=4,6)
                     ELSEIF (NVAL .EQ. 5) THEN
                        WRITE (LOUT, 3070) (PFAL(L,NFAL),L=4,8)
                     ELSE
                        WRITE (LOUT, 2020) SUB(N)(:ILEN)
                        KERR = .TRUE.
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
C
         ELSEIF (UPCASE(KEY, 3) .EQ. 'REV') THEN
C
C        REVERSE ARRHENIUS PARAMETERS
C
            IF (LFAL .OR. LREV .OR. NSPEC(II).LT.0) THEN
               KERR = .TRUE.
               IF (LFAL) WRITE (LOUT, 1090) SUB(N)(:ILEN)
               IF (LREV) WRITE (LOUT, 2050) SUB(N)(:ILEN)
               IF (NSPEC(II) .LT. 0) WRITE (LOUT, 2060) SUB(N)(:ILEN)
            ELSE
               LREV = .TRUE.
               NREV = NREV+1
               IREV(NREV) = II
               CALL IPPARR (RSTR,1,NPAR,RPAR(1,NREV),NVAL,IER,LOUT)
               IF (IER .NE. 0) KERR = .TRUE.
               WRITE (LOUT, 1900) '   Reverse Arrhenius coefficients:',
     1                           (RPAR(L,NREV),L=1,3)
            ENDIF
C
         ELSEIF (UPCASE(KEY, 3) .EQ. 'RLT') THEN
C
C        REVERSE LANDAU-TELLER PARAMETERS
C
            IF (LFAL .OR. LRLT .OR. NSPEC(II).LT.0) THEN
               KERR = .TRUE.
               IF (LFAL) WRITE (LOUT, 1070) SUB(N)(:ILEN)
               IF (LRLT) WRITE (LOUT, 2080) SUB(N)(:ILEN)
               IF (NSPEC(II) .LT. 0) WRITE (LOUT, 1080) SUB(N)(:ILEN)
            ELSE
               LRLT = .TRUE.
               NRLT = NRLT + 1
               IRLT(NRLT) = II
               CALL IPPARR (RSTR,1,NLAR,RLAN(1,NRLT),NVAL,IER,LOUT)
               IF (IER .NE. 0) KERR = .TRUE.
               WRITE (LOUT, 3040) (RLAN(L,NRLT),L=1,2)
            ENDIF
C
         ELSEIF (UPCASE(KEY, 2) .EQ. 'HV') THEN
C
C        RADIATION WAVELENGTH ENHANCEMENT FACTOR
C
            IF (.NOT.LWL) THEN
               WRITE (LOUT, 1000) SUB(N)(:ILEN)
               KERR = .TRUE.
            ELSE
               CALL IPPARR (RSTR,1,1,VAL,NVAL,IER,LOUT)
               IF (IER .EQ. 0) THEN
                  WL(NWL) = WL(NWL)*VAL(1)
                  WRITE (LOUT, 3020) ABS(WL(NWL))
               ELSE
                  WRITE (LOUT, 1000) SUB(N)(:ILEN)
                  KERR = .TRUE.
               ENDIF
            ENDIF
C
         ELSEIF (UPCASE(KEY, 2) .EQ. 'LT') THEN
C
C        LANDAU-TELLER PARAMETERS
C
            IF (LFAL .OR. LLAN) THEN
               KERR = .TRUE.
               IF (LFAL) WRITE (LOUT, 1060) SUB(N)(:ILEN)
               IF (LLAN) WRITE (LOUT, 2070) SUB(N)(:ILEN)
            ELSE
               LLAN = .TRUE.
               NLAN = NLAN + 1
               ILAN(NLAN) = II
               CALL IPPARR (RSTR,1,NLAR,PLAN(1,NLAN),NVAL,IER,LOUT)
               IF (IER .NE. 0) THEN
                  WRITE (LOUT, 1010) SUB(N)(:ILEN)
                  KERR = .TRUE.
               ENDIF
               WRITE (LOUT, 3000) (PLAN(L,NLAN),L=1,2)
            ENDIF
C
         ELSEIF (UPCASE(KEY,4).EQ.'FORD' .OR.
     1           UPCASE(KEY,4).EQ.'RORD') THEN
             LFORD = (UPCASE(KEY,4) .EQ. 'FORD')
             LRORD = (UPCASE(KEY,4) .EQ. 'RORD')
             IF (NORD.EQ.0 .OR.(NORD.GT.0 .AND. IORD(NORD).NE.II)) THEN
                NORD = NORD + 1
                IORD(NORD) = II
                NKORD = 0
C
                IF (NRNU.GT.0 .AND. IRNU(NRNU).EQ.II) THEN
                   DO 111 L = 1, 6
                      IF (NUNK(L,II) .NE. 0) THEN
                         NKORD = NKORD + 1
                         IF (RNU(L,NRNU) .LT. 0.0) THEN
                            KORD(NKORD,NORD) = -NUNK(L,II)
                            RORD(NKORD,NORD) = ABS(RNU(L,NRNU))
                         ELSE
                            KORD(NKORD,NORD) = NUNK(L,II)
                            RORD(NKORD,NORD) = RNU(L,NRNU)
                         ENDIF
                      ENDIF
  111              CONTINUE
               ELSE
                   DO 113 L = 1, 6
                      IF (NUNK(L,II) .NE. 0) THEN
                         NKORD = NKORD + 1
                         IF (NU(L,II) .LT. 0) THEN
                            KORD(NKORD,NORD) = -NUNK(L,II)
                            RORD(NKORD,NORD) =  IABS(NU(L,II))
                         ELSE
                            KORD(NKORD,NORD) = NUNK(L,II)
                            RORD(NKORD,NORD) = NU(L,II)
                         ENDIF
                      ENDIF
  113              CONTINUE
                ENDIF
             ENDIF
C
             CALL IPNPAR (RSTR, 1, ISTR, ISTART)
             IF (ISTART .GE. 1) THEN
                CALL IPPARR (ISTR, 1, 1, VAL, NVAL, IER, LOUT)
                CALL CKCOMP (RSTR(:ISTART-1), KNAME, KK, K)
                IF (LFORD) K = -K
                NK = 0
                DO 121 L = 1, MAXORD
C
                   IF (KORD(L,NORD).EQ.0) THEN
                      NK = L
                      GO TO 122
                   ELSEIF (KORD(L,NORD).EQ.K) THEN
                      IF (LFORD) THEN
                         WRITE (LOUT,*)
     1'       Warning...changing order for reactant...',
     2                   KNAME(-K)
                      ELSE
                         WRITE (LOUT,*)
     1'       Warning...changing order for product...',
     2                   KNAME(K)
                      ENDIF
                      NK = L
                      GO TO 122
                   ENDIF
  121           CONTINUE
  122           CONTINUE
                KORD(NK,NORD) = K
                RORD(NK,NORD) = VAL(1)
                IF (LFORD) THEN
                   WRITE (LOUT, 3015) KNAME(-K),VAL(1)
                ELSE
                   WRITE (LOUT, 3016) KNAME(K),VAL(1)
                ENDIF
            ENDIF
C

         ELSE
C
C        ENHANCED THIRD BODIES
C
            CALL CKCOMP (KEY, KNAME, KK, K)
            IF (K .EQ. 0) THEN
               WRITE (LOUT, 1040) KEY(:ILASCH(KEY))
               KERR = .TRUE.
            ELSE
               IF (.NOT.LTHB) THEN
                  KERR = .TRUE.
                  WRITE (LOUT, 1020) SUB(N)(:ILEN)
               ELSE
                  IF (NTBS(NTHB) .EQ. MAXTB) THEN
                     KERR = .TRUE.
                     WRITE (LOUT, 1030) SUB(N)(:ILEN)
                  ELSE
                     CALL IPPARR (RSTR, 1, 1, VAL, NVAL, IER, LOUT)
                     IF (IER .EQ. 0) THEN
                        WRITE (LOUT, 3010) KNAME(K),VAL(1)
                        NTBS(NTHB) = NTBS(NTHB) + 1
                        NKTB(NTBS(NTHB),NTHB) = K
                        AIK(NTBS(NTHB),NTHB) = VAL(1)
                     ELSE
                        WRITE (LOUT, 1020) SUB(N)(:ILEN)
                        KERR = .TRUE.
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
  500 CONTINUE
C
C     FORMATS
C
 1000 FORMAT (6X,'Error in HV declaration...',A)
 1010 FORMAT (6X,'Error in LT declaration..',A)
 1020 FORMAT (6X,'Error in third body declaration...',A)
 1030 FORMAT (6X,'Error...more than MAXTB third bodies...',A)
 1040 FORMAT (6X,'Error...undeclared species...',A)
 1050 FORMAT (6X,'Error...this is not a fall-off reaction...',A)
 1060 FORMAT (6X,'Error...LT declared in fall-off reaction...',A)
 1070 FORMAT (6X,'Error...RLT declared in fall-off reaction...',A)
 1080 FORMAT (6X,'Error...RLT declared in irreversible reaction...',A)
 1090 FORMAT (6X,'Error...REV declared in fall-off reaction...',A)
 2000 FORMAT (6X,'Error...LOW declared more than once...',A)
 2010 FORMAT (6X,'Error...TROE declared more than once...',A)
 2020 FORMAT (6X,'Error in fall-off parameters...',A)
 2030 FORMAT (6X,'Error...cannot use both TROE and SRI...',A)
 2040 FORMAT (6X,'Error...SRI declared more than once...',A)
 2050 FORMAT (6X,'Error...REV declared more than once...',A)
 2060 FORMAT (6X,'Error...REV declared for irreversible reaction...',A)
 2070 FORMAT (6X,'Error...LT declared more than once...',A)
 2080 FORMAT (6X,'Error...RLT declared more than once...',A)
 2090 FORMAT (6X,'Error in auxiliary data...',A)
 3000 FORMAT (9X,'Landau-Teller parameters: B=',E12.5,', C=',E12.5)
 3010 FORMAT (9X,A16,' Enhanced by ',1PE12.3)
 3015 FORMAT (7X,A16,' Forward order ',1PE12.3)
 3016 FORMAT (7X,A16,' Reverse order ',1PE12.3)
 3020 FORMAT (9X,'Radiation wavelength (A): ',F10.2)
C 1900 FORMAT (6X,A,T51,E10.3,F7.3,F11.3)
 1900 FORMAT (6X, A, T53, 1PE8.2, 2X, 0PF5.1, 2X, F9.1)
 3040 FORMAT (9X,'Reverse Landau-Teller parameters: B=',E12.5,
     1           ', C=',E12.5)
 3050 FORMAT (6X,'Low pressure limit:',3E13.5)
 3060 FORMAT (6X,'SRI centering:     ',3E13.5)
 3070 FORMAT (6X,'SRI centering:     ',5E13.5)
 3080 FORMAT (6X,'TROE centering:    ',3E13.5)
 3090 FORMAT (6X,'TROE centering:    ',4E13.5)
 4000 FORMAT (6X,'Declared duplicate reaction...')
      END
C----------------------------------------------------------------------C
      SUBROUTINE CKPRNT (MDIM, MAXTP, MM, ENAME, KK, KNAME, WTM,
     1                   KPHSE, KCHRG, NT, T, TLO, TMID, THI, KNCF,
     2                   ITHRM, LOUT, KERR)
C
C     Prints species interpreter output and checks for completeness.
C----------------------------------------------------------------------C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION WTM(*), KPHSE(*), KCHRG(*), T(MAXTP,*),
     1          NT(*), KNCF(MDIM,*), IPLUS(10)
      LOGICAL KERR, ITHRM(*)
      CHARACTER ENAME(*)*(*), KNAME(*)*(*), IPHSE(3)*1, INUM(10)*1
      DATA IPHSE/'S','G','L'/
      DATA INUM/'0','1','2','3','4','5','6','7','8','9'/
C
      WRITE (LOUT, 400) (ENAME(M), M = 1, MM)
      WRITE (LOUT, 300)
C
      DO 100 K = 1, KK
C
         IF (T(1,K) .LT. 0.0) T(1,K) = TLO
         IF (T(2,K) .LT. 0.0) T(2,K) = TMID
         IF (T(3,K) .LT. 0.0) T(NT(K),K) = THI
         WRITE (LOUT, 500) K, KNAME(K), IPHSE(KPHSE(K)+2), KCHRG(K),
     1                    WTM(K), T(1,K), T(NT(K),K), (KNCF(M,K),M=1,MM)
         IF (T(1,K) .GE. T(NT(K),K)) THEN
            KERR = .TRUE.
            WRITE (LOUT, 240)
         ENDIF
         IF (T(1,K) .GT. T(2,K)) THEN
            WRITE (LOUT, 250)
            KERR = .TRUE.
         ENDIF
         IF (T(NT(K),K) .LT. T(2,K)) THEN
            WRITE (LOUT, 260)
            KERR = .TRUE.
         ENDIF
C
C        each species must have thermodynamic data
C
         IF (.NOT. ITHRM(K)) THEN
            KERR = .TRUE.
            WRITE (LOUT, 200)
         ENDIF
C
C        a species cannot start with a number
C
         CALL CKCOMP (KNAME(K)(:1), INUM, 10, I)
         IF (I .GT. 0) THEN
            KERR = .TRUE.
            WRITE (LOUT, 210)
         ENDIF
C
C        if '+' sign is used in a species name,
C           examples of legal species symbols with + are:
C           OH(+)2, OH(+2), OH+, OH++, OH+++, OH(+), OH(++),
C           OH[+OH], OH2+, OH+2
C
C           examples of illegal species symbols with + are:
C           +OH        (symbol starts with a +, this will cause
C                       confusion in a reaction)
C           OH(+OH)    (symbol in parentheses is another species-
C                       this arrangement is reserved for a fall-off
C                       reaction)
C           OH+OH      (plus delimits other species names, this
C                       will cause confusion in a reaction)
C
         NPLUS = 0
         DO 50 N = 1, ILASCH(KNAME(K))
            IF (KNAME(K)(N:N) .EQ. '+') THEN
               NPLUS = NPLUS + 1
               IPLUS(NPLUS) = N
            ENDIF
   50    CONTINUE
         DO 60 N = 1, NPLUS
            I1 = IPLUS(N)
            IF (I1 .EQ. 1) THEN
               WRITE (LOUT, 220)
               KERR = .TRUE.
            ELSE
C
C              is there another species name in parentheses
C
               IF (KNAME(K)(I1-1:I1-1) .EQ. '(') THEN
                  I1 = I1 + 1
                  I2 = I1 + INDEX(KNAME(K)(I1:),')')-1
                  IF (I2 .GT. I1) THEN
                     CALL CKCOMP (KNAME(K)(I1:I2-1), KNAME, KK, KNUM)
                     IF (KNUM .GT. 0) THEN
                        WRITE (LOUT, 230)
                        KERR = .TRUE.
                     ENDIF
                  ENDIF
               ENDIF
C
C              is there another species name after a +
C
               I1 = I1 + 1
               IF (N .LT. NPLUS) THEN
                  DO 55 L = N+1, NPLUS
                     I2 = IPLUS(L)
                     IF (I2 .GT. I1) THEN
                        CALL CKCOMP (KNAME(K)(I1:I2-1),KNAME,KK,KNUM)
                        IF (KNUM .GT. 0) THEN
                           WRITE (LOUT, 230)
                           KERR = .TRUE.
                        ENDIF
                     ENDIF
   55             CONTINUE
               ENDIF
C
               I2 = ILASCH(KNAME(K))
               IF (I2 .GE. I1) THEN
                  CALL CKCOMP (KNAME(K)(I1:I2), KNAME, KK, KNUM)
                  IF (KNUM .GT. 0) THEN
                     WRITE (LOUT, 230)
                     KERR = .TRUE.
                  ENDIF
               ENDIF
            ENDIF
   60    CONTINUE
C
  100 CONTINUE
      WRITE (LOUT, 300)
      RETURN
C
  200 FORMAT (6X,'Error...no thermodynamic properties for species')
  210 FORMAT (6X,'Error...species starts with a number')
  220 FORMAT (6X,'Error...species starts with a plus')
  230 FORMAT (6X,'Error...illegal + in species name')
  240 FORMAT (6X,'Error...High temperature must be < Low temperature')
  250 FORMAT (6X,'Error...Low temperature must be <= Mid temperature')
  260 FORMAT (6X,'Error...High temperature must be => Mid temperature')
  300 FORMAT (1X,79('-'))
C  400 FORMAT (1X,79('-'),/21X,'C',/18X,'P',2X,'H',/18X,'H',2X,'A',
C     1        /18X,'A',2X,'R',/1X,'SPECIES',10X,'S',2X,'G',2X,
C     2        'MOLECULAR',3X,'TEMPERATURE',4X,'ELEMENT COUNT',/1X,
C     3        'CONSIDERED',7X,'E',2X,'E',2X,'WEIGHT',6X,'LOW',5X,
C     4        'HIGH',3X,15(A3),/1X,79('-'))
C  500 FORMAT (I4,'. ',A10,2X,A1,I3,F11.5,2(F8.1),15(I3))
C
  400 FORMAT (1X,79('-'),/T26,'C',/T24,'P H',/T24,'H A',/T24,'A R',
     1       /1X,'SPECIES',T24,'S G',T28,'MOLECULAR',T38,'TEMPERATURE',
     2       T52,'ELEMENT COUNT',
     3       /1X,'CONSIDERED',T24,'E E',T28,'WEIGHT',T38,'LOW',
     4       T45,'HIGH',T52,15(A3))
  500 FORMAT (1X,I3,'. ',A16,T24,A1,T26,I1,T28,F9.5,T38,F6.1,T45,F6.1,
     1       T51,15(I3))
      END
C----------------------------------------------------------------------C
      SUBROUTINE CPREAC (II, MAXSP, NSPEC, NPAR, PAR, RPAR, AUNITS,
     1                   EUNITS, NREAC, NUNK, NU, KCHRG, MDIM, MM, KNCF,
     2                   IDUP, NFAL, IFAL, KFAL, NFAR, PFAL, IFOP, NREV,
     3                   IREV, NTHB, ITHB, NLAN, ILAN, NRLT, IRLT, KERR,
     4                   LOUT, NRNU, IRNU, RNU, CKMIN)
C
C     Prints reaction interpreter output and checks for reaction
C     balance, duplication, and missing data in 'REV' reactions;
C     correct units of Arrhenius parameters
C
C     Input: II     - the index number of the reaction
C            MAXSP  - maximum number of species allowed in a reaction
C            NSPEC  - array of the number of species in the reactions
C            NPAR   - the number of Arrhenius parameters required
C            PAR    - matrix of Arrhenius parameters for the reactions
C            RPAR   - matrix of reverse Arrhenius parameters for the
C                     reactions which declared them
C            AUNITS - character string which describes the input units
C                     of A, the pre-exponential factor PAR(1,I)
C            EUNITS - character string which describes the input units
C                     of E, the activation energy PAR(3,I)
C            NREAC  - array of the number of reactants in the reactions
C            NUNK   - matrix of the species numbers of the reactants
C                     and products in the reactions
C            NU     - matrix of the stoichiometric coefficients of the
C                     reactants and products in the reactions
C            KCHRG  - array of the electronic charges of the species
C            MDIM   - the maximum number of elements allowed
C            MM     - the actual number of elements declared
C            KNCF   - matrix of elemental composition of the species
C            IDUP   - array of integer flags to indicate duplicate
C                     reactions
C            NFAL   - total number of reactions with fall-off
C            IFAL   - array of the NFAL reaction numbers
C            NFAR   - maximum number of fall-off parameters allowed
C            PFAL   - matrix of fall-off parameters for the NFAL
C                     reactions
C            IFOP   - array of integer fall-off types for the NFAL
C                     reactions
C            NREV   - total number of reactions with reverse parameters
C            IREV   - array of the NREV reaction numbers
C            NTHB   - total number of reactions with third-bodies
C            ITHB   - array of the NTHB reaction numbers
C            NLAN   - total number of reactions with Landauer-Teller
C                     parameters
C            ILAN   - array of the NLAN reaction numbers
C            NRLT   - total number of reactions with reverse
C                     Landauer-Teller parameters
C            IRLT   - array of the NRLT reaction numbers
C            KERR   - logical error flag
C            LOUT   - unit number for output messages
C
C----------------------------------------------------------------------C
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION NSPEC(*), PAR(NPAR,*), RPAR(NPAR,*), NREAC(*),
     1          NUNK(MAXSP,*), NU(MAXSP,*), KCHRG(*), KNCF(MDIM,*),
     2          IDUP(*), IFAL(*), KFAL(*), PFAL(NFAR,*), IFOP(*),
     3          IREV(*), ITHB(*), ILAN(*), IRLT(*), IRNU(*),
     4          RNU(MAXSP,*)
      CHARACTER*(*) AUNITS, EUNITS
      LOGICAL IERR,KERR,LREV,LLAN,LRLT
C
      IF (NRNU.GT.0 .AND. (II.EQ.IRNU(NRNU))) THEN
         CALL CKRBAL (MAXSP, NUNK(1,II), RNU(1,NRNU), MDIM, MM, KCHRG, 
     1                KNCF, CKMIN, IERR)
      ELSE
         CALL CKBAL (MAXSP, NUNK(1,II), NU(1,II), MDIM, MM, KCHRG, KNCF,
     1               IERR)
      ENDIF
C
      IF (IERR) THEN
         KERR = .TRUE.
         WRITE (LOUT, 1060)
      ENDIF
C
      CALL CKDUP (II, MAXSP, NSPEC, NREAC, NU, NUNK, NFAL, IFAL, KFAL,
     1            ISAME)
C
      IF (ISAME .GT. 0) THEN
         IF (IDUP(ISAME).NE.0 .AND. IDUP(II).NE.0) THEN
            IDUP(ISAME) = ABS(IDUP(ISAME))
            IDUP(II)    = ABS(IDUP(II))
         ELSE
            N1 = 0
            N2 = 0
            IF (NTHB .GT. 1) THEN
               DO 150 N = 1, NTHB
                  IF (ITHB(N) .EQ. ISAME) N1 = 1
                  IF (ITHB(N) .EQ. II)    N2 = 1
  150          CONTINUE
            ENDIF
            IF (N1 .EQ. N2) THEN
               KERR = .TRUE.
               WRITE (LOUT, 1050) ISAME
            ENDIF
         ENDIF
      ENDIF
C
      IF (NFAL.GT.0 .AND. IFAL(NFAL).EQ.II .AND. IFOP(NFAL).LT.0) THEN
         KERR = .TRUE.
         WRITE (LOUT, 1020)
      ENDIF
C
      LREV = (NREV.GT.0 .AND. IREV(NREV).EQ.II)
      LLAN = (NLAN.GT.0 .AND. ILAN(NLAN).EQ.II)
      LRLT = (NRLT.GT.0 .AND. IRLT(NRLT).EQ.II)
      IF (LREV .AND. LLAN .AND. (.NOT.LRLT)) THEN
         KERR = .TRUE.
         WRITE (LOUT, 1030)
      ENDIF
      IF (LRLT .AND. (.NOT.LLAN)) THEN
         KERR = .TRUE.
         WRITE (LOUT, 1040)
      ENDIF
      IF (LRLT .AND. (.NOT.LREV)) THEN
         KERR = .TRUE.
         WRITE (LOUT, 1045)
      ENDIF
C
      IF (EUNITS .EQ. 'KELV') THEN
         EFAC = 1.0
      ELSEIF (EUNITS .EQ. 'CAL/') THEN
C        convert E from cal/mole to Kelvin
         EFAC = 1.0 / 1.987
      ELSEIF (EUNITS .EQ. 'KCAL') THEN
C        convert E from kcal/mole to Kelvin
         EFAC = 1000.0 / 1.987
      ELSEIF (EUNITS .EQ. 'JOUL') THEN
C        convert E from Joules/mole to Kelvin
         EFAC = 1.0 / 8.314
      ELSEIF (EUNITS .EQ. 'KJOU') THEN
C        convert E from Kjoules/mole to Kelvin
         EFAC = 1000.0 / 8.314
      ENDIF
      PAR(3,II) = PAR(3,II) * EFAC
C
C      IF (NREV.GT.0 .AND. IREV(NREV).EQ.II) RPAR(3,II)=RPAR(3,II)*EFAC
C      IF (NFAL.GT.0 .AND. IFAL(NFAL).EQ.II) PFAL(3,II)=PFAL(3,II)*EFAC
C
      IF (NREV.GT.0 .AND. IREV(NREV).EQ.II)
     1    RPAR(3,NREV) = RPAR(3,NREV) * EFAC
      IF (NFAL.GT.0 .AND. IFAL(NFAL).EQ.II)
     1    PFAL(3,NFAL) = PFAL(3,NFAL) * EFAC
C
      IF (AUNITS .EQ. 'MOLC') THEN
         NSTOR = 0
         NSTOP = 0
         DO 50 N = 1, MAXSP
            IF (NU(N,II) .LT. 0) THEN
C              sum of stoichiometric coefficients of reactants
               NSTOR = NSTOR + ABS(NU(N,II))
            ELSEIF (NU(N,II) .GT. 0) THEN
C              sum of stoichiometric coefficients of products
               NSTOP = NSTOP + NU(N,II)
            ENDIF
   50    CONTINUE
C
         AVAG = 6.023E23
C
         IF (NFAL.GT.0 .AND. IFAL(NFAL).EQ.II) THEN
C
C           fall-off reaction, "(+M)" or "(+species name)" does not
C           count except in "LOW" A-factor;
C           reverse-rate declarations are not allowed
C
            IF (NSTOR.GT.0) PAR(1,II) = PAR(1,II) * AVAG**(NSTOR-1)
            NSTOR = NSTOR + 1
            IF (NSTOR.GT.0) PFAL(1,NFAL) = PFAL(1,NFAL)*AVAG**(NSTOR-1)
C
         ELSEIF (NTHB.GT.0 .AND. ITHB(NTHB).EQ.II) THEN
C
C           third body reaction, "+M" counts as species in
C           forward and reverse A-factor conversion
C
            NSTOR = NSTOR + 1
            NSTOP = NSTOP + 1
            IF (NSTOR.GT.0) PAR(1,II) = PAR(1,II) * AVAG**(NSTOR-1)
            IF (NREV.GT.0 .AND. IREV(NREV).EQ.II .AND. NSTOP.GT.0)
     1          RPAR(1,NREV) = RPAR(1,NREV) * AVAG**(NSTOP-1)
C
         ELSE
C
C           not third-body or fall-off reaction, but may have
C           reverse rates.
C
            IF (NSTOR .GT. 0) PAR(1,II) = PAR(1,II) * AVAG**(NSTOR-1)
            IF (NREV.GT.0 .AND. IREV(NREV).EQ.II .AND. NSTOP.GT.0)
     1          RPAR(1,NREV) = RPAR(1,NREV) * AVAG**(NSTOP-1)
         ENDIF
      ENDIF
C
 1020 FORMAT (6X,'Error...no LOW parameters given for fall-off...')
 1030 FORMAT (6X,'Error...reverse T-L required...')
 1040 FORMAT (6X,'Error...forward T-L required...')
 1045 FORMAT (6X,'Error...REV parameters must be given with RTL...')
 1050 FORMAT (6X,'Error...undeclared duplicate to reaction number ',I3)
 1060 FORMAT (6X,'Error...reaction does not balance...')
      RETURN
      END
C----------------------------------------------------------------------C
      SUBROUTINE CKBAL (MXSPEC, KSPEC, KCOEF, MDIM, MM, KCHRG, KNCF,
     1                  IERR)
C
C     Checks elemental balance of reactants vs. products.
C     Checks charge balance of reaction.
C
C     Input:  MXSPEC - number of species allowed in a reaction
C             KSPEC(N),N=1,MXSPEC- array of species numbers in reaction
C             KCOEF(N) - stoichiometric coefficients of the species
C             MDIM  - maximum number of elements allowed
C             MM    - actual integer number of elements
C             KCHRG(K) - ionic charge Kth species
C             KNCF(M,K)- integer elemental composition of Kth species
C     Output: KERR  - logical, =.TRUE. if reaction does not balance
C                                      F. Rupley, Div. 8245, 5/13/86
C----------------------------------------------------------------------C
C*****precision > double
       IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION KSPEC(*), KCOEF(*), KNCF(MDIM,*), KCHRG(*)
      LOGICAL IERR
C
      IERR = .FALSE.
C
C     charge balance
C
      KBAL = 0
      DO 50 N = 1, MXSPEC
         IF (KSPEC(N) .NE. 0) 
     1   KBAL = KBAL + KCOEF(N)*KCHRG(KSPEC(N))
   50 CONTINUE
      IF (KBAL .NE. 0) IERR = .TRUE.
C
C     element balance
C
      DO 100 M = 1, MM
         MBAL = 0
         DO 80 N = 1, MXSPEC
            IF (KSPEC(N) .NE. 0)
     1      MBAL = MBAL + KCOEF(N)*KNCF(M,KSPEC(N))
   80    CONTINUE
         IF (MBAL .NE. 0) IERR = .TRUE.
  100 CONTINUE
      RETURN
      END
C----------------------------------------------------------------------C
      SUBROUTINE CKRBAL (MXSPEC, KSPEC, RCOEF, MDIM, MM, KCHRG, KNCF,
     1                   CKMIN, IERR)
C
C     Checks elemental balance of reactants vs. products.
C     Checks charge balance of reaction.
C
C     Input:  MXSPEC - number of species allowed in a reaction
C             KSPEC(N),N=1,MXSPEC- array of species numbers in reaction
C             RCOEF(N) - stoichiometric coefficients of the species
C             MDIM  - maximum number of elements allowed
C             MM    - actual integer number of elements
C             KCHRG(K) - ionic charge Kth species
C             KNCF(M,K)- integer elemental composition of Kth species
C     Output: KERR  - logical, =.TRUE. if reaction does not balance
C                                      F. Rupley, Div. 8245, 5/13/86
C----------------------------------------------------------------------C
C*****precision > double
       IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION KSPEC(*), RCOEF(*), KNCF(MDIM,*), KCHRG(*)
      LOGICAL IERR
C
      IERR = .FALSE.
C
C     charge balance
C
      SBAL = 0
      DO 50 N = 1, MXSPEC
         IF (KSPEC(N) .NE. 0)
     1   SBAL = SBAL + RCOEF(N)*KCHRG(KSPEC(N))
   50 CONTINUE
      IF (ABS(SBAL) .GT. CKMIN) IERR = .TRUE.
C
C     element balance
C
      DO 100 M = 1, MM
         SMBAL = 0
         DO 80 N = 1, MXSPEC
            IF (KSPEC(N) .NE. 0)
     1      SMBAL = SMBAL + RCOEF(N)*KNCF(M,KSPEC(N))
   80    CONTINUE
         IF (ABS(SMBAL) .GT. CKMIN) IERR = .TRUE.
  100 CONTINUE
      RETURN
      END
C----------------------------------------------------------------------C
      SUBROUTINE CKDUP (I, MAXSP, NS, NR, NU, NUNK, NFAL, IFAL, KFAL,
     1                  ISAME)
C
C     Checks reaction I against the (I-1) reactions for duplication
C----------------------------------------------------------------------C
C*****precision > double
       IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      DIMENSION NS(*), NR(*), NU(MAXSP,*), NUNK(MAXSP,*), IFAL(*),
     1          KFAL(*)
C
      ISAME = 0
      NRI = NR(I)
      NPI = ABS(NS(I)) - NR(I)
C
      DO 500 J = 1, I-1
C
         NRJ = NR(J)
         NPJ = ABS(NS(J)) - NR(J)
C
         IF (NRJ.EQ.NRI .AND. NPJ.EQ.NPI) THEN
C
            NSAME = 0
            DO 20 N = 1, MAXSP
               KI = NUNK(N,I)
               NI = NU(N,I)
C
               DO 15 L = 1, MAXSP
                  KJ = NUNK(L,J)
                  NJ = NU(L,J)
                  IF (NJ.NE.0 .AND. KJ.EQ.KI .AND. NJ.EQ.NI)
     1            NSAME = NSAME + 1
   15          CONTINUE
   20       CONTINUE
C
            IF (NSAME .EQ. ABS(NS(J))) THEN
C
C           same products, reactants, coefficients, check fall-off
C           third body
C
               IF (NFAL.GT.0 .AND. IFAL(NFAL).EQ.I) THEN
                  DO 22 N = 1, NFAL-1
                     IF (J.EQ.IFAL(N) .AND. KFAL(N).EQ.KFAL(NFAL)) THEN
                        ISAME = J
                        RETURN
                     ENDIF
   22             CONTINUE
                  RETURN
               ENDIF
C
               ISAME = J
               RETURN
            ENDIF
         ENDIF
C
         IF (NPI.EQ.NRJ .AND. NPJ.EQ.NRI) THEN
C
            NSAME = 0
            DO 30 N = 1, MAXSP
               KI = NUNK(N,I)
               NI = NU(N,I)
C
               DO 25 L = 1, MAXSP
                  KJ = NUNK(L,J)
                  NJ = NU(L,J)
                  IF (NJ.NE.0 .AND. KJ.EQ.KI .AND. -NJ.EQ.NI)
     1            NSAME = NSAME + 1
   25          CONTINUE
   30       CONTINUE
C
            IF (NSAME.EQ.ABS(NS(J)) .AND.
     1          (NS(J).GT.0 .OR. NS(I).GT.0)) THEN
C
C           same products as J reactants, and vice-versa
C
               IF (NFAL.GT.0 .AND. IFAL(NFAL).EQ.I) THEN
                  DO 32 N = 1, NFAL-1
                     IF (J.EQ.IFAL(N) .AND. KFAL(N).EQ.KFAL(NFAL)) THEN
                        ISAME = J
                        RETURN
                     ENDIF
   32             CONTINUE
                  RETURN
               ENDIF
C
               ISAME = J
               RETURN
            ENDIF
         ENDIF
C
  500 CONTINUE
      RETURN
      END
C----------------------------------------------------------------------C
      SUBROUTINE CKISUB (LINE, SUB, NSUB)
C
C     Generates an array of CHAR*(*) substrings from a CHAR*(*) string,
C     using blanks or tabs as delimiters
C
C     Input:  LINE  - a CHAR*(*) line
C     Output: SUB   - a CHAR*(*) array of substrings
C             NSUB  - number of substrings found
C     A '!' will comment out a line, or remainder of the line.
C                                      F. Rupley, Div. 8245, 5/15/86
C----------------------------------------------------------------------C
C*****precision > double
       IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      CHARACTER*(*) SUB(*), LINE
      NSUB = 0
C
      DO 5 N = 1, LEN(LINE)
        IF (ICHAR(LINE(N:N)) .EQ. 9) LINE(N:N) = ' '
    5 CONTINUE
C
      IF (IPPLEN(LINE) .LE. 0) RETURN
C
      ILEN = ILASCH(LINE)
C
      NSTART = IFIRCH(LINE)
   10 CONTINUE
      ISTART = NSTART
      NSUB = NSUB + 1
      SUB(NSUB) = ' '
C
      DO 100 I = ISTART, ILEN
         ILAST = INDEX(LINE(ISTART:),' ') - 1
         IF (ILAST .GT. 0) THEN
            ILAST = ISTART + ILAST - 1
         ELSE
            ILAST = ILEN
         ENDIF
         SUB(NSUB) = LINE(ISTART:ILAST)
         IF (ILAST .EQ. ILEN) RETURN
C
         NSTART = ILAST + IFIRCH(LINE(ILAST+1:))
C
C        Does SUB have any slashes?
C
         I1 = INDEX(SUB(NSUB),'/')
         IF (I1 .LE. 0) THEN
            IF (LINE(NSTART:NSTART) .NE. '/') GO TO 10
            NEND = NSTART + INDEX(LINE(NSTART+1:),'/')
            IND = INDEX(SUB(NSUB),' ')
            SUB(NSUB)(IND:) = LINE(NSTART:NEND)
            IF (NEND .EQ. ILEN) RETURN
            NSTART = NEND + IFIRCH(LINE(NEND+1:))
            GO TO 10
         ENDIF
C
C        Does SUB have 2 slashes?
C
         I2 = INDEX(SUB(NSUB)(I1+1:),'/')
         IF (I2 .GT. 0) GO TO 10
C
         NEND = NSTART + INDEX(LINE(NSTART+1:),'/')
         IND = INDEX(SUB(NSUB),' ') + 1
         SUB(NSUB)(IND:) = LINE(NSTART:NEND)
         IF (NEND .EQ. ILEN) RETURN
         NSTART = NEND + IFIRCH(LINE(NEND+1:))
C        GO TO 10
  100 CONTINUE
      RETURN
      END
C----------------------------------------------------------------------C
      SUBROUTINE IPNPAR (LINE, NPAR, IPAR, ISTART)
C
C     Returns CHAR*(*) IPAR substring of CHAR*(*) string LINE which
C     contains NPAR real parameters
C
C     Input:     LINE - a CHAR*(*) line
C                NPAR - number of parameters expected
C     Output:    IPAR - the substring of parameters only
C                ISTART - the starting location of IPAR substring
C     A '!' will comment out a line, or remainder of the line.
C                                      F. Rupley, Div. 8245, 5/14/86
C----------------------------------------------------------------------C
C*****precision > double
       IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      CHARACTER*(*) LINE,IPAR
C
C----------Find Comment String (! signifies comment)
C
      ILEN = IPPLEN(LINE)
      ISTART = 0
      N = 0
      IF (ILEN.GT.0) THEN
         DO 40 I = ILEN, 1, -1
            ISTART = I
            IPAR = ' '
            IPAR = LINE(ISTART:ILEN)
            IF (LINE(I:I).NE.' ') THEN
               IF (I .EQ. 1) RETURN
               IF (LINE(I-1:I-1) .EQ. ' ') THEN
                  N = N + 1
                  IF (N .EQ. NPAR) RETURN
               ENDIF
            ENDIF
   40    CONTINUE
      ENDIF
      RETURN
      END
C----------------------------------------------------------------------C
      SUBROUTINE IPPARI(STRING, ICARD, NEXPEC, IVAL, NFOUND, IERR, LOUT)
C   BEGIN PROLOGUE  IPPARI
C   REFER TO  IPGETI
C   DATE WRITTEN  850625   (YYMMDD)
C   REVISION DATE 851725   (YYMMDD)
C   CATEGORY NO.  J3.,J4.,M2.
C   KEYWORDS  PARSE
C   AUTHOR  CLARK,G.L.,GROUP C-3 LOS ALAMOS NAT'L LAB
C   PURPOSE  Parses integer variables from a character variable.  Called
C            by IPGETI, the IOPAK routine used for interactive input.
C   DESCRIPTION
C
C-----------------------------------------------------------------------
C  IPPARI may be used for parsing an input record that contains integer
C  values, but was read into a character variable instead of directly
C  into integer variables.
C  The following benefits are gained by this approach:
C    - specification of only certain elements of the array is allowed,
C      thus letting the others retain default values
C    - variable numbers of values may be input in a record, up to a
C      specified maximum
C    - control remains with the calling program in case of an input
C      error
C    - diagnostics may be printed by IPPARI to indicate the nature
C      of input errors
C
C   The contents of STRING on input indicate which elements of IVAL
C   are to be changed from their entry values, and values to which
C   they should be changed on exit.  Commas and blanks serve as
C   delimiters, but multiple blanks are treated as a single delimeter.
C   Thus, an input record such as:
C     '   1,   2,,40000   , ,60'
C   is interpreted as the following set of instructions by IPGETR:
C
C     (1) set IVAL(1) = 1
C     (2) set IVAL(2) = 2
C     (3) leave IVAL(3) unchanged
C     (4) set IVAL(4) = 40000
C     (5) leave IVAL(5) unchanged
C     (6) set IVAL(6) = 60
C
C   IPPARI will print diagnostics on the default output device, if
C   desired.
C
C   IPPARI is part of IOPAK, and is written in ANSI FORTRAN 77
C
C   Examples:
C
C      Assume IVAL = (0, 0, 0) and NEXPEC = 3 on entry:
C
C   input string           IVAL on exit            IERR    NFOUND
C   -------------          ----------------------  ----    ------
C  '  2 ,   3 45 '         (2, 3, 45)                0       3
C  '2.15,,3'               (2, 0, 3)                 1       0
C  '3X, 25, 2'             (0, 0, 0)                 1       0
C  '10000'                 (10000, 0, 0)             2       1
C
C      Assume IVAL = (0, 0, 0, 0) and NEXPEC = -4 on entry:
C
C   input string           IVAL on exit            IERR    NFOUND
C   -------------          ----------------------  ----    ------
C  '1, 2'                  (1, 2)                    0       2
C  ',,37  400'             (0, 0, 37, 400)           0       4
C  ' 1,,-3,,5'             (1, 0, -3, 0)             3       4
C
C  arguments: (I=input,O=output)
C  -----------------------------
C  STRING (I) - the character string to be parsed.
C
C  ICARD  (I) - data statement number, and error processing flag
C         < 0 : no error messages printed
C         = 0 : print error messages, but not ICARD
C         > 0 : print error messages, and ICARD
C
C  NEXPEC (I) - number of real variables expected to be input.  If
C         < 0, the number is unknown, and any number of values
C         between 0 and abs(nexpec) may be input.  (see NFOUND)
C
C  PROMPT (I) - prompting string, character type.  A question
C         mark will be added to form the prompt at the screen.
C
C  IVAL (I,O) - the integer value or values to be modified.  On entry,
C       the values are printed as defaults.  The formal parameter
C       corresponding to IVAL must be dimensioned at least NEXPEC
C       in the calling program if NEXPEC > 1.
C
C  NFOUND (O) - the number of real values represented in STRING,
C         only in the case that there were as many or less than
C         NEXPEC.
C
C  IERR (O) - error flag:
C       = 0 if no errors found
C       = 1 syntax errors or illegal values found
C       = 2 for too few values found (NFOUND < NEXPEC)
C       = 3 for too many values found (NFOUND > NEXPEC)
C-----------------------------------------------------------------------
C
C   REFERENCES  (NONE)
C   ROUTINES CALLED  IFIRCH,ILASCH
C   END PROLOGUE  IPPARI
C*****precision > double
       IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
C
      CHARACTER STRING*(*), ITEMP*80
      DIMENSION IVAL(*)
      CHARACTER *8 FMT(14)
      LOGICAL OKINCR
C
C   FIRST EXECUTABLE STATEMENT  IPPARI
      IERR   = 0
      NFOUND = 0
      NEXP = IABS(NEXPEC)
      IE = ILASCH(STRING)
      IF (IE .EQ. 0) GO TO 500
      NC = 1
C
C--- OKINCR is a flag that indicates it's OK to increment
C--- NFOUND, the index of the array into which the value
C--- should be read.  It is set false when a space follows
C--- an integer value substring, to keep incrementing from
C--- occurring if a comma should be encountered before the
C--- next value.
C
      OKINCR = .TRUE.
C
C--- begin overall loop on characters in string
C
100   CONTINUE
C
      IF (STRING(NC:NC) .EQ. ',') THEN
         IF (OKINCR .OR. NC .EQ. IE) THEN
            NFOUND = NFOUND + 1
         ELSE
            OKINCR = .TRUE.
         ENDIF
C
         GO TO 450
      ENDIF
      IF (STRING(NC:NC) .EQ. ' ') GO TO 450
C
C--- first good character (non-delimeter) found - now find
C--- last good character
C
      IBS = NC
160   CONTINUE
      NC = NC + 1
      IF (NC .GT. IE) GO TO 180
      IF (STRING(NC:NC) .EQ. ' ')THEN
         OKINCR = .FALSE.
      ELSEIF (STRING(NC:NC) .EQ. ',')THEN
         OKINCR = .TRUE.
      ELSE
         GO TO 160
      ENDIF
C
C--- end of substring found - read value into integer array
C
180   CONTINUE
      NFOUND = NFOUND + 1
      IF (NFOUND .GT. NEXP) THEN
         IERR = 3
         GO TO 500
      ENDIF
C
      IES = NC - 1
      NCH = IES - IBS + 1
      DATA FMT/' (I1)', ' (I2)', ' (I3)', ' (I4)', ' (I5)',
     1   ' (I6)', ' (I7)', ' (I8)', ' (I9)', '(I10)',
     2   '(I11)', '(I12)', '(I13)', '(I14)'/
      ITEMP = ' '
      ITEMP = STRING(IBS:IES)
      READ (ITEMP(1:NCH), FMT(NCH), ERR = 400) IVAL(NFOUND)
      GO TO 450
400   CONTINUE
      IERR = 1
      GO TO 510
450   CONTINUE
      NC = NC + 1
      IF (NC .LE. IE) GO TO 100
C
500   CONTINUE
      IF (NEXPEC .GT. 0 .AND. NFOUND .LT. NEXP) IERR = 2
510   CONTINUE
C
      IF (IERR .EQ. 0 .OR. ICARD .LT. 0)RETURN
      IF (ICARD .NE. 0) WRITE (LOUT, '(A,I3)')
     1   '!! ERROR IN DATA STATEMENT NUMBER', ICARD
      IF (IERR .EQ. 1)
     1    WRITE (LOUT, '(A)')'SYNTAX ERROR, OR ILLEGAL VALUE'
      IF (IERR .EQ. 2) WRITE (LOUT, '(A,I2, A, I2)')
     1   ' TOO FEW DATA ITEMS.  NUMBER FOUND = ' , NFOUND,
     2   '  NUMBER EXPECTED = ', NEXPEC
      IF (IERR .EQ. 3) WRITE (LOUT, '(A,I2)')
     1   ' TOO MANY DATA ITEMS.  NUMBER EXPECTED = ', NEXPEC
      END
C
      SUBROUTINE IPPARR(STRING, ICARD, NEXPEC, RVAL, NFOUND, IERR, LOUT)
C   BEGIN PROLOGUE  IPPARR
C   REFER TO  IPGETR
C   DATE WRITTEN  850625   (YYMMDD)
C   REVISION DATE 851625   (YYMMDD)
C   CATEGORY NO.  J3.,J4.,M2.
C   KEYWORDS  PARSE
C   AUTHOR  CLARK,G.L.,GROUP C-3 LOS ALAMOS NAT'L LAB
C   PURPOSE  Parses real variables from a character variable.  Called
C            by IPGETR, the IOPAK routine used for interactive input.
C   DESCRIPTION
C
C-----------------------------------------------------------------------
C  IPPARR may be used for parsing an input record that contains real
C  values, but was read into a character variable instead of directly
C  into real variables.
C  The following benefits are gained by this approach:
C    - specification of only certain elements of the array is allowed,
C      thus letting the others retain default values
C    - variable numbers of values may be input in a record, up to a
C      specified maximum
C    - control remains with the calling program in case of an input
C      error
C    - diagnostics may be printed by IPPARR to indicate the nature
C      of input errors
C
C   The contents of STRING on input indicate which elements of RVAL
C   are to be changed from their entry values, and values to which
C   they should be changed on exit.  Commas and blanks serve as
C   delimiters, but multiple blanks are treated as a single delimeter.
C   Thus, an input record such as:
C     '   1.,   2,,4.e-5   , ,6.e-6'
C   is interpreted as the following set of instructions by IPGETR:
C
C     (1) set RVAL(1) = 1.0
C     (2) set RVAL(2) = 2.0
C     (3) leave RVAL(3) unchanged
C     (4) set RVAL(4) = 4.0E-05
C     (5) leave RVAL(5) unchanged
C     (6) set RVAL(6) = 6.0E-06
C
C   IPPARR will print diagnostics on the default output device, if
C   desired.
C
C   IPPARR is part of IOPAK, and is written in ANSI FORTRAN 77
C
C   Examples:
C
C      Assume RVAL = (0., 0., 0.) and NEXPEC = 3 on entry:
C
C   input string           RVAL on exit            IERR    NFOUND
C   -------------          ----------------------  ----    ------
C  '  2.34e-3,  3 45.1'    (2.34E-03, 3.0, 45.1)     0       3
C  '2,,3.-5'               (2.0, 0.0, 3.0E-05)       0       3
C  ',1.4,0.028E4'          (0.0, 1.4, 280.0)         0       3
C  '1.0, 2.a4, 3.0'        (1.0, 0.0, 0.0)           1       1
C  '1.0'                   (1.0, 0.0, 0.0)           2       1
C
C      Assume RVAL = (0.,0.,0.,0.) and NEXPEC = -4 on entry:
C
C   input string           RVAL on exit            IERR    NFOUND
C   -------------          ----------------------  ----    ------
C  '1.,2.'                 (1.0, 2.0)                0       2
C  ',,3  4.0'              (0.0, 0.0, 3.0, 4.0)      0       4
C  '1,,3,,5.0'             (0.0, 0.0, 3.0, 0.0)      3       4
C
C  arguments: (I=input,O=output)
C  -----------------------------
C  STRING (I) - the character string to be parsed.
C
C  ICARD  (I) - data statement number, and error processing flag
C         < 0 : no error messages printed
C         = 0 : print error messages, but not ICARD
C         > 0 : print error messages, and ICARD
C
C  NEXPEC (I) - number of real variables expected to be input.  If
C         < 0, the number is unknown, and any number of values
C         between 0 and abs(nexpec) may be input.  (see NFOUND)
C
C  PROMPT (I) - prompting string, character type.  A question
C         mark will be added to form the prompt at the screen.
C
C  RVAL (I,O) - the real value or values to be modified.  On entry,
C       the values are printed as defaults.  The formal parameter
C       corresponding to RVAL must be dimensioned at least NEXPEC
C       in the calling program if NEXPEC > 1.
C
C  NFOUND (O) - the number of real values represented in STRING,
C         only in the case that there were as many or less than
C         NEXPEC.
C
C  IERR (O) - error flag:
C       = 0 if no errors found
C       = 1 syntax errors or illegal values found
C       = 2 for too few values found (NFOUND < NEXPEC)
C       = 3 for too many values found (NFOUND > NEXPEC)
C-----------------------------------------------------------------------
C
C   REFERENCES  (NONE)
C   ROUTINES CALLED  IFIRCH,ILASCH
C   END PROLOGUE  IPPARR
C*****precision > double
       IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      CHARACTER STRING*(*), ITEMP*80
      DIMENSION RVAL(*)
      CHARACTER *8 FMT(22)
      LOGICAL OKINCR
C
C   FIRST EXECUTABLE STATEMENT  IPPARR
      IERR   = 0
      NFOUND = 0
      NEXP = IABS(NEXPEC)
      IE = ILASCH(STRING)
      IF (IE .EQ. 0) GO TO 500
      NC = 1
C
C--- OKINCR is a flag that indicates it's OK to increment
C--- NFOUND, the index of the array into which the value
C--- should be read.  It is set negative when a space follows
C--- a real value substring, to keep incrementing from
C--- occurring if a comma should be encountered before the
C--- next value.
C
      OKINCR = .TRUE.
C
C--- begin overall loop on characters in string
C
100   CONTINUE
C
      IF (STRING(NC:NC) .EQ. ',') THEN
         IF (OKINCR) THEN
            NFOUND = NFOUND + 1
         ELSE
            OKINCR = .TRUE.
         ENDIF
C
         GO TO 450
      ENDIF
      IF (STRING(NC:NC) .EQ. ' ') GO TO 450
C
C--- first good character (non-delimeter) found - now find
C--- last good character
C
      IBS = NC
160   CONTINUE
      NC = NC + 1
      IF (NC .GT. IE) GO TO 180
      IF (STRING(NC:NC) .EQ. ' ')THEN
         OKINCR = .FALSE.
      ELSEIF (STRING(NC:NC) .EQ. ',')THEN
         OKINCR = .TRUE.
      ELSE
         GO TO 160
      ENDIF
C
C--- end of substring found - read value into real array
C
180   CONTINUE
      NFOUND = NFOUND + 1
      IF (NFOUND .GT. NEXP) THEN
         IERR = 3
         GO TO 500
      ENDIF
C
      DATA FMT/     ' (E1.0)', ' (E2.0)', ' (E3.0)', ' (E4.0)',
     1   ' (E5.0)', ' (E6.0)', ' (E7.0)', ' (E8.0)', ' (E9.0)',
     2   '(E10.0)', '(E11.0)', '(E12.0)', '(E13.0)', '(E14.0)',
     3   '(E15.0)', '(E16.0)', '(E17.0)', '(E18.0)', '(E19.0)',
     4   '(E20.0)', '(E21.0)', '(E22.0)'/
      IES = NC - 1
      NCH = IES - IBS + 1
      ITEMP = ' '
      ITEMP = STRING(IBS:IES)
      READ (ITEMP(1:NCH), FMT(NCH), ERR = 400) RVAL(NFOUND)
      GO TO 450
400   CONTINUE
      WRITE (LOUT, 555) STRING(IBS:IES)
  555 FORMAT (A)
      IERR = 1
      GO TO 510
450   CONTINUE
      NC = NC + 1
      IF (NC .LE. IE) GO TO 100
C
500   CONTINUE
      IF (NEXPEC .GT. 0 .AND. NFOUND .LT. NEXP) IERR = 2
510   CONTINUE
C
      IF (IERR .EQ. 0 .OR. ICARD .LT. 0) RETURN
      IF (ICARD .NE. 0) WRITE (LOUT, '(A,I3)')
     1   '!! ERROR IN DATA STATEMENT NUMBER', ICARD
      IF (IERR .EQ. 1)
     1   WRITE (LOUT, '(A)')'SYNTAX ERROR, OR ILLEGAL VALUE'
      IF (IERR .EQ. 2) WRITE (LOUT, '(A,I2, A, I2)')
     1   ' TOO FEW DATA ITEMS.  NUMBER FOUND = ' , NFOUND,
     2   '  NUMBER EXPECTED = ', NEXPEC
      IF (IERR .EQ. 3) WRITE (LOUT, '(A,I2)')
     1   ' TOO MANY DATA ITEMS.  NUMBER EXPECTED = ', NEXPEC
      END
C
      FUNCTION IFIRCH(STRING)
C   BEGIN PROLOGUE  IFIRCH
C   DATE WRITTEN   850626
C   REVISION DATE  850626
C   CATEGORY NO.  M4.
C   KEYWORDS  CHARACTER STRINGS,SIGNIFICANT CHARACTERS
C   AUTHOR  CLARK,G.L.,GROUP C-3 LOS ALAMOS NAT'L LAB
C   PURPOSE  Determines first significant (non-blank) character
C            in character variable
C   DESCRIPTION
C
C-----------------------------------------------------------------------
C  IFIRCH locates the first non-blank character in a string of
C  arbitrary length.  If no characters are found, IFIRCH is set = 0.
C  When used with the companion routine ILASCH, the length of a string
C  can be determined, and/or a concatenated substring containing the
C  significant characters produced.
C-----------------------------------------------------------------------
C
C   REFERENCES  (NONE)
C   ROUTINES CALLED  (NONE)
C   END PROLOGUE IFIRCH
C*****precision > double
       IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      CHARACTER* (*)STRING
C
C   FIRST EXECUTABLE STATEMENT IFIRCH
      NLOOP = LEN(STRING)
C
      IF (NLOOP .EQ. 0) THEN
         IFIRCH = 0
         RETURN
      ENDIF
C
      DO 100 I = 1, NLOOP
         IF (STRING(I:I) .NE. ' ') GO TO 120
100   CONTINUE
C
      IFIRCH = 0
      RETURN
120   CONTINUE
      IFIRCH = I
      END
      FUNCTION ILASCH(STRING)
C   BEGIN PROLOGUE  ILASCH
C   DATE WRITTEN   850626
C   REVISION DATE  850626
C   CATEGORY NO.  M4.
C   KEYWORDS  CHARACTER STRINGS,SIGNIFICANT CHARACTERS
C   AUTHOR  CLARK,G.L.,GROUP C-3 LOS ALAMOS NAT'L LAB
C   PURPOSE  Determines last significant (non-blank) character
C            in character variable
C   DESCRIPTION
C
C-----------------------------------------------------------------------
C  IFIRCH locates the last non-blank character in a string of
C  arbitrary length.  If no characters are found, ILASCH is set = 0.
C  When used with the companion routine IFIRCH, the length of a string
C  can be determined, and/or a concatenated substring containing the
C  significant characters produced.
C  Note that the FORTRAN intrinsic function LEN returns the length
C  of a character string as declared, rather than as filled.  The
C  declared length includes leading and trailing blanks, and thus is
C  not useful in generating 'significant' substrings.
C-----------------------------------------------------------------------
C
C   REFERENCES  (NONE)
C   ROUTINES CALLED  (NONE)
C   END PROLOGUE IFIRCH
C*****precision > double
       IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      CHARACTER*(*) STRING
C
C***FIRST EXECUTABLE STATEMENT ILASCH
      NLOOP = LEN(STRING)
      IF (NLOOP.EQ.0) THEN
         ILASCH = 0
         RETURN
      ENDIF
C
      DO 100 I = NLOOP, 1, -1
         IF (STRING(I:I) .NE. ' ') GO TO 120
100   CONTINUE
C
120   CONTINUE
      ILASCH = I
      END
C----------------------------------------------------------------------C
C
      SUBROUTINE CKCOMP (IST, IRAY, II, I)
C
C  START PROLOGUE
C
C  SUBROUTINE CKCOMP (IST, IRAY, II, I)*
C     Returns the index of an element of a reference character
C     string array which corresponds to a character string;
C     leading and trailing blanks are ignored.
C
C
C  INPUT
C     IST   - A character string.
C                  Data type - CHARACTER*(*)
C     IRAY  - An array of character strings;
C             dimension IRAY(*) at least II
C                  Data type - CHARACTER*(*)
C     II    - The length of IRAY.
C                  Data type - integer scalar.
C
C  OUTPUT
C     I     - The first integer location in IRAY in which IST
C             corresponds to IRAY(I); if IST is not also an
C             entry in IRAY, I=0.
C
C  END PROLOGUE
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      CHARACTER*(*) IST, IRAY(*)
C
      I = 0
      DO 10 N = II, 1, -1
         IS1 = IFIRCH(IST)
         IS2 = ILASCH(IST)
         IR1 = IFIRCH(IRAY(N))
         IR2 = ILASCH(IRAY(N))
         IF ( IS2.GE.IS1 .AND. IS2.GT.0 .AND.
     1        IR2.GE.IR1 .AND. IR2.GT.0 .AND.
     2        IST(IS1:IS2).EQ.IRAY(N)(IR1:IR2) ) I=N
   10 CONTINUE
      RETURN
      END
C
C----------------------------------------------------------------------C
      SUBROUTINE CKUNIT (LINE, AUNITS, EUNITS, IUNITS)
C
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
      CHARACTER*(*) LINE, IUNITS, AUNITS, EUNITS
      CHARACTER*4 UPCASE
C
      AUNITS = ' '
      EUNITS = ' '
      IUNITS = ' '
      LCHAR = ILASCH(LINE)
      DO 85 N = 1, ILASCH(LINE)-3
         IND = ILASCH(IUNITS)
         IF (EUNITS .EQ. ' ') THEN
            IF (UPCASE(LINE(N:), 4)     .EQ. 'CAL/') THEN
               EUNITS = 'CAL/'
               IF (IUNITS .EQ. ' ') THEN
                  IUNITS = 'E units cal/mole'
               ELSE
                  IUNITS(IND:) = ', E units cal/mole'
               ENDIF
            ELSEIF (UPCASE(LINE(N:), 4) .EQ. 'KCAL') THEN
               EUNITS = 'KCAL'
               IF (IUNITS .EQ. ' ') THEN
                  IUNITS = 'E units Kcal/mole'
               ELSE
                  IUNITS(IND:) = ', E units Kcal/mole'
               ENDIF
            ELSEIF (UPCASE(LINE(N:), 4) .EQ. 'JOUL') THEN
               EUNITS = 'JOUL'
               IF (IUNITS .EQ. ' ') THEN
                  IUNITS = 'E units Joules/mole'
               ELSE
                  IUNITS(IND:) = ', E units Joules/mole'
               ENDIF
            ELSEIF (UPCASE(LINE(N:), 4) .EQ. 'KJOU') THEN
               EUNITS = 'KJOU'
               IF (IUNITS .EQ. ' ') THEN
                  IUNITS = 'E units Kjoule/mole'
               ELSE
                  IUNITS(IND:) = ', E units Kjoule/mole'
               ENDIF
            ELSEIF (UPCASE(LINE(N:), 4) .EQ. 'KELV') THEN
               EUNITS = 'KELV'
               IF (IUNITS .EQ. ' ') THEN
                  IUNITS = 'E units Kelvins'
               ELSE
                  IUNITS(IND:) = ', E units Kelvins'
               ENDIF
            ENDIF
         ENDIF
         IF (AUNITS .EQ. ' ') THEN
            IF (UPCASE(LINE(N:), 4) .EQ. 'MOLE') THEN
               IF (N+4.LE.ILASCH(LINE) .AND. 
     1                    UPCASE(LINE(N+4:),1).EQ.'C') THEN
C
                  AUNITS = 'MOLC'
                  IF (IUNITS .EQ. ' ') THEN
                     IUNITS = 'A units molecules'
                  ELSE
                      IUNITS(IND:) = ', A units molecules'
                  ENDIF
               ELSE
                  AUNITS = 'MOLE'
                  IF (IUNITS .EQ. ' ') THEN
                     IUNITS = 'A units mole-cm-sec-K'
                  ELSE
                     IUNITS(IND:) = ', A units mole-cm-sec-K'
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
   85 CONTINUE
C
      IF (AUNITS .EQ. ' ') THEN
         AUNITS = 'MOLE'
         IND = ILASCH(IUNITS) + 1
         IF (IND .GT. 1) THEN
            IUNITS(IND:) = ', A units mole-cm-sec-K'
         ELSE
            IUNITS(IND:) = ' A units mole-cm-sec-K'
         ENDIF
      ENDIF
C
      IF (EUNITS .EQ. ' ') THEN
         EUNITS = 'CAL/'
         IND = ILASCH(IUNITS) + 1
         IF (IND .GT. 1) THEN
            IUNITS(IND:) = ', E units cal/mole'
         ELSE
            IUNITS(IND:) = ' E units cal/mole'
         ENDIF
      ENDIF
C
      RETURN
      END
C
C----------------------------------------------------------------------C
C
      INTEGER FUNCTION IPPLEN (LINE)
C
C  BEGIN PROLOGUE
C
C  FUNCTION IPPLEN (LINE)
C     Returns the effective length of a character string, i.e.,
C     the index of the last character before an exclamation mark (!)
C     indicating a comment.
C
C  INPUT
C     LINE  - A character string.
C
C  OUTPUT
C     IPPLEN - The effective length of the character string.
C
C  END PROLOGUE
C
C*****precision > double
       IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C      IMPLICIT REAL (A-H,O-Z), INTEGER (I-N)
C*****END precision > single
C
      CHARACTER LINE*(*)
C
      IN = IFIRCH(LINE)
      IF (IN.EQ.0 .OR. LINE(IN:IN).EQ.'!') THEN
         IPPLEN = 0
      ELSE
         IN = INDEX(LINE,'!')
         IF (IN .EQ. 0) THEN
            IPPLEN = ILASCH(LINE)
         ELSE
            IPPLEN = ILASCH(LINE(:IN-1))
         ENDIF
      ENDIF
      RETURN
      END
C
      CHARACTER*(*) FUNCTION UPCASE(ISTR, ILEN)
      CHARACTER ISTR*(*), LCASE(26)*1, UCASE(26)*1
      DATA LCASE /'a','b','c','d','e','f','g','h','i','j','k','l','m',
     1            'n','o','p','q','r','s','t','u','v','w','x','y','z'/,
     2     UCASE /'A','B','C','D','E','F','G','H','I','J','K','L','M',
     3            'N','O','P','Q','R','S','T','U','V','W','X','Y','Z'/
C
      UPCASE = ' '
      UPCASE = ISTR(:ILEN)
      JJ = MIN (LEN(UPCASE), LEN(ISTR), ILEN)
      DO 10 J = 1, JJ
         DO 10 N = 1,26
            IF (ISTR(J:J) .EQ. LCASE(N)) UPCASE(J:J) = UCASE(N)
   10 CONTINUE
      RETURN
      END
