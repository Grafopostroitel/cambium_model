! http://kiwi.atmos.colostate.edu/fortran/docs/fortran2012.key-6.pdf
MODULE double
    IMPLICIT NONE
    INTEGER, PARAMETER :: DP = KIND(0.0D0)
END MODULE double

MODULE parameters
    USE double
    IMPLICIT NONE
    REAL(KIND=DP), PARAMETER :: ONE = 1.0_DP
    INTEGER, PARAMETER :: RING = 1, CAMBIUM = 2, ENLARGEMENT = 3, MATURATION = 4
    INTEGER, PARAMETER :: DISTANCE_CELL = 4, CONCENTRATION = 11, SIZE = 8, &
        INHIBITOR = 10, AGE_CELL = 1, TIME_SPEND_CAMBIUM = 2, &
        DISTANCE_DIFFERENTIATION = 3, TIME_SWITCH_DEATH_STAGE = 5, & 
        INITIAL_CELL_SIZE = 6, WALL_THICKNESS = 9
        
END MODULE parameters

MODULE cambium_model
    USE double
    USE parameters
    IMPLICIT NONE
CONTAINS
    SUBROUTINE camb70r(ch, b, number_cells, state)
    IMPLICIT NONE
        ! P27.10.2001/7.1.2016

        !  ch(:, :)  
        !   1 - the age of the cell - возраст клетки
        !   2 - the time spend in cambium
        !   3 - the distance of differentiation
        !   4 - the distance of cell http://www.biologydiscussion.com/botany/study-notes-on-coniferopsida-botany/19182
        !   5 - the time switch to death stage
        !   6 - the initial cell size 
        !   7 - the ?
        !   8 - the size
        !   9 - the wall thickness - толщина стенки
        !  10 - the inhibitor - замедлитель 
        !  11 - the concentration - концентрация
        !  12 - ?
        ! http://biology.ru/course/content/chapter9/section1/paragraph7/theory.html#.VvJjKtKLS00
        ! http://www.activestudy.info/mikroskopicheskoe-stroenie-drevesiny-serdceviny-i-kory/ 
        ! http://www.zoodrug.ru/topic1783.html Как размножаются клетки  
        REAL(KIND=DP), DIMENSION(:, :), INTENT(INOUT) :: ch ! characteristics cells
        REAL(KIND=DP), DIMENSION(:), INTENT(IN) :: b ! cambial parameters
        INTEGER, DIMENSION(:), INTENT(INOUT) :: number_cells ! number cells in the ring, cambium, enlargement, maturation
        CHARACTER(LEN=*), INTENT(INOUT) :: state(:) ! state of cell (c, e, m, d)
        
        CHARACTER(LEN=1) :: w
        INTEGER :: i, i1, j, kk, ndiv
        REAL(KIND=DP) :: cdist, t0, v, v0, x
        
        ! Calculating some values needed below  Вычисление некоторых значений необходимых ниже
        v0 = b(26) ! Growth rate during S(synthesis), G2, and M(mitosis) phases (default 1.0) ! http://biofile.ru/bio/21539.html
        t0 = 0.1_DP ! the time-step
        ! Initial conditions Начальные условия  
        DO kk = 1, 10
            ndiv = 0 ! number of divisions per day число делений в сутки
            j = 1 ! position of cell
120 &  ! Loop of pozition from label 120 to 150
            ! characteristics of each cell характеристики каждой клетки
            w = state(j)
            SELECT CASE (w)
            CASE('c') ! cambium cells
                cdist = 0.0D0
                ! the distance cell j from initial cell
                DO i = 1, j
                    cdist = cdist + ch(8, i)
                END DO
                ch(11, j) = ch(10, j) / ch(8, j) ! concentration  
                v = rate_division(b, ch, j)
                ch(4, j) = cdist
                ch(5, j) = v
                IF(ch(8,j).LE.b(10))THEN 
                    ! cell phase G1
                    IF(ch(11,j).LT.b(24))THEN         
                        ! if concentration too low - switch cell to zone of enlargement
                        state(j) = 'e'                 
                        ch(3, j) = cdist
                        ch(6, j) = ch(8, j)
                        ch(7, j) = ch(2, j)
                        CONTINUE 
                    ELSE 
                        ! growth of cambial cell in G1
                        ch(8, j) = ch(8, j) + v * t0 ! size
                        ch(2, j) = ch(2, j) + t0
                    ENDIF
                    GOTO 150 ! next cell
                ELSEIF (ch(8, j) .LE. b(13)) THEN  
                    ! growth of cambial cell in S, G2 and M 
                    ch(8, j) = ch(8, j) + v0 * t0 
                    ch(2, j) = ch(2, j) + t0
                    GOTO 150 ! next cell
                END IF
                !  Division of cambial cell in position j
                ndiv = ndiv + 1
                ! The shift cells on one pozition 
                DO i = 1000, j + 2, -1
                    state(i) = state(i - 1)
                    DO i1 = 1, 11
                        ch(i1, i) = ch(i1, i - 1)
                    END DO
                END DO
                ! first daughter cell - position j
                ch(8, j) = ch(8, j) / 2.0 ! becomes half at cell division
                x = ch(10, j)
                ch(10, j) = x * b(25)
                IF (j.EQ.1) THEN
                    ch(10, j) = b(27) ! for initial cell y = 1
                END IF
                ch(11, j) =ch(10, j) / ch(8, j)
                ch(4, j) = ch(4, j) - ch(8, j)
                ! second daughter cell - position j + 1
                j = j + 1
                ch(8, j) = ch(8, j - 1)
                ch(10, j) = x * (1 - b(25))
                ch(2, j) = ch(2, j - 1)
                IF(j.EQ.2) ch(2, j) = 0.0
                ch(4, j) = ch(4, j - 1) + ch(8, j)
                ch(11, j) = ch(10, j) / ch(8, j)          
                state(j) = 'c' ! state is still cambium
            CASE DEFAULT
            END SELECT
150 &
            j = j + 1 
            IF ((j <= (number_cells(RING) + ndiv)) .AND. (j <= 1000)) THEN 
                GOTO 120 ! processes next cell  
            END IF
            ! determination of the number cells in each zone - определение количества ячеек в каждой зоне
            number_cells(RING)        = number_cells(RING) + ndiv
            number_cells(CAMBIUM)     = 0
            number_cells(ENLARGEMENT) = 0
            number_cells(MATURATION)  = 0
            DO i = 1, number_cells(RING) 
                w = state(i)
                SELECT CASE (w)
                    CASE ('c')
                        number_cells(CAMBIUM) = number_cells(CAMBIUM) + 1 ! number in cambium
                    CASE ('e')
                        number_cells(ENLARGEMENT) = number_cells(ENLARGEMENT) + 1 ! number in enlargement 
                    CASE ('m')
                        number_cells(MATURATION) = number_cells(MATURATION) + 1 ! number in enlargement 
                    CASE DEFAULT
                END SELECT
            END DO
        END DO 
    END SUBROUTINE camb70r

    FUNCTION rate_division(b, ch, j) RESULT (vp)
        IMPLICIT NONE
        REAL(KIND=DP), DIMENSION(:), INTENT(IN) :: b
        REAL(KIND=DP), DIMENSION(:, :), INTENT(IN) :: ch
        INTEGER, INTENT(IN) :: j
        
        REAL(KIND=DP) :: vp
        
        vp = b(22) + b(21) * EXP(-b(23) * ch(CONCENTRATION, j))
        IF(vp <= 0.0_DP) THEN 
            vp = 0.0_DP ! vp cannot exceed a horizontal line
        END IF
        
    END FUNCTION rate_division
END MODULE cambium_model

! Cambium Model Run
PROGRAM growth_cambium 
USE cambium_model
USE parameters
IMPLICIT NONE

    INTEGER, PARAMETER :: PTIME = 300
    REAL(KIND=DP)      :: ch(12, 1000) ! characteristics cells
    REAL(KIND=DP)      :: outx(500), outy(500)
    ! The control of division of cambium cells (20)
    ! The cellular constants (20)
    REAL(KIND=DP)      :: b(100) 
    ! number cells in the ring, cambium, enlargement, maturation
    INTEGER            :: number_cells(4) 
    INTEGER            :: nt(PTIME), nc(PTIME), ne(PTIME)
    INTEGER            :: op(10)
    INTEGER            :: Time
    CHARACTER(LEN=1)   :: state(1000)
    CHARACTER(LEN=41)  :: legx, legy
    
    INTEGER            :: i, j, k
    REAL(KIND=DP)      :: xd, xmax, ymax, z

    ! DEBUG WRITE
    OPEN(99, FILE='dbg.out')
    ! Reading of the parameters and the options of cambium model
    OPEN(10, STATUS='OLD', FILE='growth_cam7.1.par')
    READ(10, *) 
    READ(10, *)
    READ(10, '(I4)') (op(I), i = 1, 10)
    DO j = 0, 1
        READ(10, *) 
        DO i = 1, 20
            READ(10, '(F10.3)') b(i + j * 20)
        END DO
    END DO
    CLOSE(10)

    ! The initial conditions
    ch     = 0.0_DP 
    state  = 'z'
    state(1) = 'c'
    
    number_cells(RING)        = 1
    number_cells(CAMBIUM)     = 1
    number_cells(ENLARGEMENT) = 0
    number_cells(MATURATION)  = 0

    ch(SIZE, 1)            = b(13) / 1.4_DP
    ch(WALL_THICKNESS, 1)  = b(7)
    ch(INHIBITOR, 1)       = b(27)
    ch(CONCENTRATION, 1)   = 1.0_DP / ch(8, 1)

    OPEN(11, FILE = 'out.dat')
    WRITE(11, '(2A6, 4A9)') 'Time', 'J', 'Dist',' Inhib', 'Conc', 'Size'

    ! Window for growth rate from position (win 1)
    legx = ' Dist'
    legy = ' V(z)'
    xmax = 100
    ymax = b(22) + b(21)
    WRITE(legy(8:15), '(F5.3)') ymax
        
    OPEN(21, FILE='rate_from_position.dat')
    WRITE(21,*) '# Window for growth rate from position (Окно №1)'
    WRITE(21,*) '# The distance of cell, V(z)'
    DO Time = 1, PTIME
    
        CALL camb70r(ch, b, number_cells, state)
        WRITE(99, *) 'number_cells= ', number_cells
        nt(Time) = number_cells(RING)
        nc(Time) = number_cells(CAMBIUM)
        ne(Time) = number_cells(ENLARGEMENT)

        k = nc(Time)
        !print * , k
        DO i = 1, 10
            outx(i) = ch(DISTANCE_CELL, i)
            ! ?
            outy(i) = b(22) + b(21) * EXP(-b(23) * ch(CONCENTRATION, i))
            WRITE(21, *) outx(i), outy(i)
            !WRITE(21, *) 1.0_DP, 1.0_DP
        END DO
        DO i = 1, k
            !outx(i) = ch(DISTANCE_CELL, i)
            !outy(i) = b(22) + b(21) * EXP(-b(23) * ch(CONCENTRATION, i))
            !WRITE(21, *) outx(i), outy(i)
            WRITE(11, '(2I6, 4F9.4 )') Time, i, ch(DISTANCE_CELL, i), & 
                ch(INHIBITOR, i), ch(CONCENTRATION, i), ch(SIZE, i)
        END DO 
    END DO ! Time
    CLOSE(11)
    CLOSE(21)
    CLOSE(99)

    ! Graph the dynamic cambium cells (win 2)
    xmax=300
    ymax=20
    DO i = 1, PTIME
        IF(Nc(i).GT.ymax)ymax=Nc(i)
    END DO
    legx='Day'
    legy='Nc' ! Number of cells in the cambial zone
    WRITE(legy(8:15),'(F4.0)')ymax
    
    CALL create_gnuplot_script('dynamic_cambium_cells.dat', 'dynamic_cambium_cells.plt', &
        'Graph the dynamic cambium cells (Окно №2)', 'Day', 'Nc')

    OPEN(21, FILE='dynamic_cambium_cells.dat')
    WRITE(21,*) '# Graph the dynamic cambium cells (win 2)'
    WRITE(21,*) '# day Nc'
    k = PTIME
    DO i = 1, k
        outx(i) = i
        outy(i) = nc(i)
        WRITE(21, *) outx(i), outy(i)
    ENDDO
    CLOSE(21)
    
    ! Graph the age of cells (win 3)
    xmax=300
    ymax=20
    DO i = 2, PTIME
        IF(ch(2,i).GT.ymax)ymax=ch(2,I)
    ENDDO
    legx='Cell'
    legy='Age,'
    WRITE(legy(8:15),'(F4.0)')ymax
    
    CALL create_gnuplot_script('age_of_cells.dat', 'age_of_cells.plt', &
        'Graph the age of cells (win 3)', 'Cell', 'Age')
    
    OPEN(21, FILE='age_of_cells.dat')
    WRITE(21,*) '# Graph the age of cells (win 3)'
    WRITE(21,*) '# Cell Age'
    OPEN(22, FILE='age_of_cells.plt')
    k=200
    DO i=1,k
        outx(i)=i+1
        outy(i)=ch(2,i+1)
        WRITE(21, *) outx(i), outy(i)
    ENDDO
    CLOSE(21)

    ! The graph of position of cells (win 4 )
    xmax=300
    ymax=20
    DO i = 2, PTIME
        IF(ch(4,i).GT.ymax)ymax=ch(4,I)
    ENDDO
    legx='Cell'
    legy='Dist,'
    WRITE(legy(8:15),'(F4.0)')ymax
    
    CALL create_gnuplot_script('position_of_cells.dat', 'position_of_cells.plt', &
        'The graph of position of cells (win 4 )', legx, legy)
    
    OPEN(21, FILE='position_of_cells.dat')
    WRITE(21,*) '# The graph of position of cells (win 4 )'
    WRITE(21,*) '# Cell Dist'
    OPEN(22, FILE='position_of_cells.plt')
    k=200
    DO i=1,k
        outx(i)=i+1
        outy(i)=ch(4,i+1)
        WRITE(21, *) outx(i), outy(i)
    ENDDO
    CLOSE(21)
    CLOSE(22)

    ! The graph of concentration of inhibitor in cells (win 5)
    xmax=300
    ymax=0.0
    DO i=12,200
        IF(ch(11,i).GT.ymax)ymax=ch(11,i)
    ENDDO
    legx='Cell'
    legy='Conc,'
    WRITE(legy(8:15),'(F4.0)')ymax
    
    CALL create_gnuplot_script('inhibitor_in_cells.dat', 'inhibitor_in_cells.plt', &
        'The graph of concentration of inhibitor in cells (win 5)', legx, legy)

    OPEN(21, FILE='inhibitor_in_cells.dat')
    WRITE(21,*) '# The graph of concentration of inhibitor in cells (win 5)'
    WRITE(21,*) '# Cell Conc'
    OPEN(22, FILE='inhibitor_in_cells.plt')
    k = 200
    DO i = 1, k
        outx(i) = i + 59
        outy(i) = ch(11, i + 59)
        WRITE(21, *) outx(i), outy(i)
    ENDDO
    CLOSE(21)
    CLOSE(22)

    OPEN(10, FILE='cambium.dat')
    WRITE(10, '(4A6)') 'Time','N','Nc','N-Nc'
    DO i = 1, PTIME
        WRITE(10, '(4I6)') i, nt(i), nc(i), ne(i)
    END DO
    CLOSE(10)
    
    OPEN(10, FILE='ch.dat')
    WRITE(10, '(A5,11A7)')'J',' ','age','Dist','Z',' ','Do',' ','Size',' ','Y','Conc'
    DO i = 1, PTIME
        WRITE(10, '(I5,11F7.2)') i, (ch(j, i), j = 1, 11)
    END DO
    CLOSE(10)
    
    ! The relationship the growth rate on concentration (win6)
    xmax = b(27) / 5
    ymax = (b(22) + b(21)) * 1.2
    legx='Concentration, ?, xmax= '
    WRITE(legx(8:15), '(F5.2)') xmax
    legy='V(c), ?, ymax= '
    WRITE(legy(8:15), '(F5.2)') ymax
    
    CALL create_gnuplot_script('rate_on_concentration.dat', 'rate_on_concentration.plt', &
        'The relationship the growth rate on concentration (win №6)', legx, legy)

    OPEN(21, FILE='rate_on_concentration.dat')
    WRITE(21,*) '# The relationship the growth rate on concentration ( win )'
    WRITE(21,*) '# Cell Conc'
    OPEN(22, FILE='rate_on_concentration.plt')
    k = 100
    xd = xmax / 100
    DO i = 1, k
        z = xd * i
        outx(i) = z
        outy(i) = (b(22) + b(21) * EXP(-b(23) * z))
        WRITE(21, *) outx(i), outy(i)
    ENDDO
    CLOSE(21)
    CLOSE(22)

CONTAINS

    SUBROUTINE create_gnuplot_script(xy_filename, script_filename, &
        title, legx, legy) 
    IMPLICIT NONE
        CHARACTER(LEN=*), INTENT(IN) :: xy_filename, script_filename, title, legx, legy  

        OPEN(UNIT=11, FILE=script_filename)
        
        Write(11, *) '# gnuplot cambium graph'

        Write(11, *) 'reset'

        Write(11, *) 'set terminal png size 1024,768 '
        Write(11, *) 'set output '//'"'//xy_filename//'.png"'

        Write(11, *) '#color definitions'
        Write(11, *) 'set border linewidth 1.5'
        Write(11, *) 'set style line 1 lc rgb "#0060ad" lt 1 lw 2 # --- blue'

        Write(11, *) 'unset key'

        Write(11, *) '# Axes'
        Write(11, *) 'set style line 11 lc rgb "#808080" lt 1'
        Write(11, *) 'set border 3 back ls 11'
        Write(11, *) 'set tics nomirror out scale 0.75'
        Write(11, *) '# Grid'
        Write(11, *) 'set style line 12 lc rgb "#808080" lt 0 lw 1'
        Write(11, *) 'set grid back ls 12'

        Write(11, *) 'set title '//'"'//title//'"'
        Write(11, *) 'set xlabel '//'"'//legx//'"'
        Write(11, *) 'set ylabel '//'"'//legy//'"'

        Write(11, *) 'plot  '//'"'//xy_filename//'"'//' w l ls 1'
        
        CLOSE(11)
    END SUBROUTINE create_gnuplot_script
END PROGRAM growth_cambium




! FORTRAN LINKS
! http://www.math.spbu.ru/user/rus/cluster/Doc/Library/fortran95/langref/langr_oglav.shtml
! http://ssd.sscc.ru/old/school/2008s/files/fortran.pdf
! http://www.phy.ohiou.edu/~elster/phys5071/extras/F95_primer.pdf
! http://www.cs.mtu.edu/~shene/COURSES/cs201/NOTES/F90-Subprograms.pdf
! http://www.fortran90.org/src/best-practices.html
! GnuPlot Links
! http://lowrank.net/gnuplot/
!
! c:\22\MinGW\bin\g95.exe -Wall -c "%f"
! Cambium Cell Links
! http://botweb.uwsp.edu/anatomy/
! http://dendro.cnre.vt.edu/forestbiology/cambium2_no_scene_1.swf
! http://e-lib.gasu.ru/eposobia/papina/bolprak/R_4_2.html Тема: Анатомическое строение стебля 
! https://books.google.ru/books?id=4q8oIla3fBIC&pg=PA442&lpg=PA442&dq=cambium+minimum+cell+size&source=bl&ots=BtdT3ljprv&sig=qW11OsgZVMzY6vZl7Vp7Yf_VbPk&hl=en&sa=X&redir_esc=y#v=onepage&q=cambium%20minimum%20cell%20size&f=false
! https://books.google.ru/books?id=MjvaCQAAQBAJ&pg=PT10&lpg=PT10&dq=%D1%80%D0%BE%D1%81%D1%82+%D0%BA%D0%B0%D0%BC%D0%B1%D0%B8%D1%8F+%D1%81%D0%BE%D1%81%D0%BD%D1%8B&source=bl&ots=FxHQM6Wg9c&sig=2YycWpKgT6H-TQrzN1tKVpD6U6g&hl=en&sa=X&redir_esc=y#v=onepage&q=%D1%80%D0%BE%D1%81%D1%82%20%D0%BA%D0%B0%D0%BC%D0%B1%D0%B8%D1%8F%20%D1%81%D0%BE%D1%81%D0%BD%D1%8B&f=false
! http://www.uic.edu/classes/bios/bios100/lectf03am/treetrunk.jpg
! гр. kytos — клетка

!http://www.wsl.ch/info/mitarbeitende/fonti/XCELL/Prietzkof_Dendrochronologia_2014.pdf http://www.wsl.ch/info/mitarbeitende/fonti/XCELL/index_EN http://www.wsl.ch/fe/landschaftsdynamik/dendroecology/Publikationen/Cuny_etal_New_Phytologist_2014

!! Cell Growth
! http://www.studmed.ru/docs/document3439?view=2
! http://dok.opredelim.com/docs/index-56414.html
! http://meduniver.com/Medical/genetika/fazi_kletochnogo_cikla.html
