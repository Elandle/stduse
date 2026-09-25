submodule(stduse) stduse_character
    implicit none

    contains

            !> \brief Reads the next line from a file.
        !!
        !! Reminder: in Fortran, files are streams.
        !! So if a line is read from a file (eg, by using the `read` intrinsic)
        !! that line is consumed in the stream: the next time a line is read,
        !! it is actually the next line in the file (to go back, an intrinsic
        !! such as `backspace` or `rewind` may be used).
        !!
        !! Reads the next line from file with unit `funit`, storing it in `line`
        !! with length `maxlen`. Leading spaces are removed, and `line` is padded
        !! with trailing spaces to fill up to `maxlen` length if necessary.
        !!
        !! Example: if the line being read is: </p>
        !! `nstab = 4` </p>
        !! Then if `readln` is called with `maxlen = 12`, the read in line is `"nstab = 4    "`. 
        !!
        !! \param[in]    funit  (`integer`)                      Unit of file to read from.
        !! \param[out]   iostat (`integer`)                      Status of the `read` statement on the file read from.
        !! \param[inout] line   (`character(len=:), allocatable`) String to allocate (`line` must be allocatable) read line into.
        !! \param[in]    maxlen (`integer, optional`)            Maximum length of a line to read and the exit length of `line`. Lines with data to be read in longer than `maxlen` being read into may result in errors. Default value: `100`.
        module procedure readln
            ! Cannot edit maxlen since it is optional, so using another integer to hold the value used
            integer :: actual_maxlen

            if (present(maxlen)) then
                actual_maxlen = maxlen
            else
                ! default maximum length of 1024
                actual_maxlen = 1024
            endif

            ! allocate if not already allocated
            if (.not. allocated(line)) then
                allocate(character(actual_maxlen) :: line)
            endif

            ! force line to be maxlen long
            if (len(line) .ne. actual_maxlen) then
                deallocate(line)
                allocate(character(actual_maxlen) :: line)
            endif

            read(unit=funit, fmt="(a)", iostat=iostat) line
            line = ltrim(line)
        endprocedure readln

        !> \brief Removes leading and trailing blank spaces from an input string.
        !!
        !! Example: if `line = "  abc de   "`, then `ltrim(line)` returns
        !! `"abc de"`.
        !!
        !! \param[in] line   (`character(len=*)`)       String to remove leading and trailing blank spaced from.
        !! \result    ltrim  (`character(len=trimlen)`) `line` with leading and trailing blank spaces removed, with a length `trimlen` of `line` minus the amount of leading and trailing blank spaces.
        module procedure ltrim
            ltrim = trim(adjustl(line))
        endprocedure ltrim
endsubmodule stduse_character