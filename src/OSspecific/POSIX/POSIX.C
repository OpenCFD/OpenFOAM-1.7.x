/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    POSIX versions of the functions declared in OSspecific.H

\*---------------------------------------------------------------------------*/

#ifdef solarisGcc
# define _SYS_VNODE_H
#endif

#include "OSspecific.H"
#include "POSIX.H"
#include "foamVersion.H"
#include "fileName.H"
#include "fileStat.H"
#include "timer.H"

#include <fstream>
#include <cstdlib>
#include <cctype>

#include <stdio.h>
#include <unistd.h>
#include <dirent.h>
#include <pwd.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/socket.h>
#include <netdb.h>

#include <netinet/in.h>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::POSIX, 0);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

pid_t Foam::pid()
{
    return getpid();
}

pid_t Foam::ppid()
{
    return getppid();
}

pid_t Foam::pgid()
{
    return getpgrp();
}

bool Foam::env(const word& envName)
{
    return getenv(envName.c_str()) != NULL;
}


Foam::string Foam::getEnv(const word& envName)
{
    char* env = getenv(envName.c_str());

    if (env)
    {
        return string(env);
    }
    else
    {
        // Return null-constructed string rather than string::null
        // to avoid cyclic dependencies in the construction of globals
        return string();
    }
}


bool Foam::setEnv
(
    const word& envName,
    const string& value,
    const bool overwrite
)
{
    return setenv(envName.c_str(), value.c_str(), overwrite) == 0;
}


Foam::word Foam::hostName()
{
    char buffer[256];
    gethostname(buffer, 256);

    return buffer;
}


Foam::word Foam::userName()
{
    struct passwd* pw = getpwuid(getuid());

    if (pw != NULL)
    {
        return pw->pw_name;
    }
    else
    {
        return word::null;
    }
}


// use $HOME environment variable or passwd info
Foam::fileName Foam::home()
{
    char* env = getenv("HOME");

    if (env != NULL)
    {
        return fileName(env);
    }
    else
    {
        struct passwd* pw = getpwuid(getuid());

        if (pw != NULL)
        {
            return pw->pw_dir;
        }
        else
        {
            return fileName::null;
        }
    }
}


Foam::fileName Foam::home(const word& userName)
{
    struct passwd* pw;

    if (userName.size())
    {
        pw = getpwnam(userName.c_str());
    }
    else
    {
        char* env = getenv("HOME");

        if (env != NULL)
        {
            return fileName(env);
        }

        pw = getpwuid(getuid());
    }

    if (pw != NULL)
    {
        return pw->pw_dir;
    }
    else
    {
        return fileName::null;
    }
}


Foam::fileName Foam::cwd()
{
    char buf[255];
    if (getcwd(buf, 255))
    {
        return buf;
    }
    else
    {
        FatalErrorIn("Foam::cwd()")
            << "Couldn't get the current working directory"
            << exit(FatalError);

        return fileName::null;
    }
}


bool Foam::chDir(const fileName& dir)
{
    return chdir(dir.c_str()) != 0;
}


Foam::fileName Foam::findEtcFile(const fileName& name, bool mandatory)
{
    // Search user files:
    // ~~~~~~~~~~~~~~~~~~
    fileName searchDir = home()/".OpenFOAM";
    if (isDir(searchDir))
    {
        // Check for user file in ~/.OpenFOAM/VERSION
        fileName fullName = searchDir/FOAMversion/name;
        if (isFile(fullName))
        {
            return fullName;
        }

        // Check for version-independent user file in ~/.OpenFOAM
        fullName = searchDir/name;
        if (isFile(fullName))
        {
            return fullName;
        }
    }


    // Search site files:
    // ~~~~~~~~~~~~~~~~~~
    searchDir = getEnv("WM_PROJECT_INST_DIR");
    if (isDir(searchDir))
    {
        // Check for site file in $WM_PROJECT_INST_DIR/site/VERSION
        fileName fullName = searchDir/"site"/FOAMversion/name;
        if (isFile(fullName))
        {
            return fullName;
        }

        // Check for version-independent site file in $WM_PROJECT_INST_DIR/site
        fullName = searchDir/"site"/name;
        if (isFile(fullName))
        {
            return fullName;
        }
    }

    // Search installation files:
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~
    searchDir = getEnv("WM_PROJECT_DIR");
    if (isDir(searchDir))
    {
        // Check for shipped OpenFOAM file in $WM_PROJECT_DIR/etc
        fileName fullName = searchDir/"etc"/name;
        if (isFile(fullName))
        {
            return fullName;
        }
    }

    // Not found
    // abort if the file is mandatory, otherwise return null
    if (mandatory)
    {
        cerr<< "--> FOAM FATAL ERROR in Foam::findEtcFile() :"
               " could not find mandatory file\n    '"
            << name.c_str() << "'\n\n" << std::endl;
        ::exit(1);
    }

    // Return null-constructed fileName rather than fileName::null
    // to avoid cyclic dependencies in the construction of globals
    return fileName();
}


bool Foam::mkDir(const fileName& pathName, mode_t mode)
{
    // empty names are meaningless
    if (pathName.empty())
    {
        return false;
    }

    // Construct instance path directory if does not exist
    if (::mkdir(pathName.c_str(), mode) == 0)
    {
        // Directory made OK so return true
        return true;
    }
    else
    {
        switch (errno)
        {
            case EPERM:
            {
                FatalErrorIn("Foam::mkDir(const fileName&, mode_t)")
                    << "The filesystem containing " << pathName
                    << " does not support the creation of directories."
                    << exit(FatalError);

                return false;
            }

            case EEXIST:
            {
                // Directory already exists so simply return true
                return true;
            }

            case EFAULT:
            {
                FatalErrorIn("Foam::mkDir(const fileName&, mode_t)")
                    << "" << pathName
                    << " points outside your accessible address space."
                    << exit(FatalError);

                return false;
            }

            case EACCES:
            {
                FatalErrorIn("Foam::mkDir(const fileName&, mode_t)")
                    << "The parent directory does not allow write "
                       "permission to the process,"<< nl
                    << "or one of the directories in " << pathName
                    << " did not allow search (execute) permission."
                    << exit(FatalError);

                return false;
            }

            case ENAMETOOLONG:
            {
                FatalErrorIn("Foam::mkDir(const fileName&, mode_t)")
                    << "" << pathName << " is too long."
                    << exit(FatalError);

                return false;
            }

            case ENOENT:
            {
                // Part of the path does not exist so try to create it
                if (pathName.path().size() && mkDir(pathName.path(), mode))
                {
                    return mkDir(pathName, mode);
                }
                else
                {
                    FatalErrorIn("Foam::mkDir(const fileName&, mode_t)")
                        << "Couldn't create directory " << pathName
                        << exit(FatalError);

                    return false;
                }
            }

            case ENOTDIR:
            {
                FatalErrorIn("Foam::mkDir(const fileName&, mode_t)")
                    << "A component used as a directory in " << pathName
                    << " is not, in fact, a directory."
                    << exit(FatalError);

                return false;
            }

            case ENOMEM:
            {
                FatalErrorIn("Foam::mkDir(const fileName&, mode_t)")
                    << "Insufficient kernel memory was available to make "
                       "directory " << pathName << '.'
                    << exit(FatalError);

                return false;
            }

            case EROFS:
            {
                FatalErrorIn("Foam::mkDir(const fileName&, mode_t)")
                    << "" << pathName
                    << " refers to a file on a read-only filesystem."
                    << exit(FatalError);

                return false;
            }

            case ELOOP:
            {
                FatalErrorIn("Foam::mkDir(const fileName&, mode_t)")
                    << "Too many symbolic links were encountered in resolving "
                    << pathName << '.'
                    << exit(FatalError);

                return false;
            }

            case ENOSPC:
            {
                FatalErrorIn("Foam::mkDir(const fileName&, mode_t)")
                    << "The device containing " << pathName
                    << " has no room for the new directory or "
                    << "the user's disk quota is exhausted."
                    << exit(FatalError);

                return false;
            }

            default:
            {
                FatalErrorIn("Foam::mkDir(const fileName&, mode_t)")
                    << "Couldn't create directory " << pathName
                    << exit(FatalError);

                return false;
            }
        }
    }
}


// Set the file mode
bool Foam::chMod(const fileName& name, const mode_t m)
{
    return ::chmod(name.c_str(), m) == 0;
}


// Return the file mode
mode_t Foam::mode(const fileName& name)
{
    fileStat fileStatus(name);
    if (fileStatus.isValid())
    {
        return fileStatus.status().st_mode;
    }
    else
    {
        return 0;
    }
}


// Return the file type: FILE or DIRECTORY
Foam::fileName::Type Foam::type(const fileName& name)
{
    mode_t m = mode(name);

    if (S_ISREG(m))
    {
        return fileName::FILE;
    }
    else if (S_ISDIR(m))
    {
        return fileName::DIRECTORY;
    }
    else
    {
        return fileName::UNDEFINED;
    }
}


// Does the name exist in the filing system?
bool Foam::exists(const fileName& name, const bool checkGzip)
{
    return mode(name) || isFile(name, checkGzip);
}


// Does the directory exist?
bool Foam::isDir(const fileName& name)
{
    return S_ISDIR(mode(name));
}


// Does the file exist?
bool Foam::isFile(const fileName& name, const bool checkGzip)
{
    return S_ISREG(mode(name)) || (checkGzip && S_ISREG(mode(name + ".gz")));
}


// Return size of file
off_t Foam::fileSize(const fileName& name)
{
    fileStat fileStatus(name);
    if (fileStatus.isValid())
    {
        return fileStatus.status().st_size;
    }
    else
    {
        return -1;
    }
}


// Return time of last file modification
time_t Foam::lastModified(const fileName& name)
{
    fileStat fileStatus(name);
    if (fileStatus.isValid())
    {
        return fileStatus.status().st_mtime;
    }
    else
    {
        return 0;
    }
}


// Read a directory and return the entries as a string list
Foam::fileNameList Foam::readDir
(
    const fileName& directory,
    const fileName::Type type,
    const bool filtergz
)
{
    // Initial filename list size
    // also used as increment if initial size found to be insufficient
    static const int maxNnames = 100;

    if (POSIX::debug)
    {
        Info<< "readDir(const fileName&, const fileType, const bool filtergz)"
            << " : reading directory " << directory << endl;
    }

    // Setup empty string list MAXTVALUES long
    fileNameList dirEntries(maxNnames);

    // Pointers to the directory entries
    DIR *source;
    struct dirent *list;

    // Temporary variables and counters
    label nEntries = 0;

    // Attempt to open directory and set the structure pointer
    if ((source = opendir(directory.c_str())) == NULL)
    {
        dirEntries.setSize(0);

        if (POSIX::debug)
        {
            Info<< "readDir(const fileName&, const fileType, "
                   "const bool filtergz) : cannot open directory "
                << directory << endl;
        }
    }
    else
    {
        // Read and parse all the entries in the directory
        while ((list = readdir(source)) != NULL)
        {
            fileName fName(list->d_name);

            // ignore files begining with ., i.e. '.', '..' and '.*'
            if (fName.size() && fName[0] != '.')
            {
                word fExt = fName.ext();

                if
                (
                    (type == fileName::DIRECTORY)
                 ||
                    (
                        type == fileName::FILE
                     && fName[fName.size()-1] != '~'
                     && fExt != "bak"
                     && fExt != "BAK"
                     && fExt != "old"
                     && fExt != "save"
                    )
                )
                {
                    if ((directory/fName).type() == type)
                    {
                        if (nEntries >= dirEntries.size())
                        {
                            dirEntries.setSize(dirEntries.size() + maxNnames);
                        }

                        if (filtergz && fExt == "gz")
                        {
                            dirEntries[nEntries++] = fName.lessExt();
                        }
                        else
                        {
                            dirEntries[nEntries++] = fName;
                        }
                    }
                }
            }
        }

        // Reset the length of the entries list
        dirEntries.setSize(nEntries);

        closedir(source);
    }

    return dirEntries;
}


// Copy, recursively if necessary, the source to the destination
bool Foam::cp(const fileName& src, const fileName& dest)
{
    // Make sure source exists.
    if (!exists(src))
    {
        return false;
    }

    fileName destFile(dest);

    // Check type of source file.
    if (src.type() == fileName::FILE)
    {
        // If dest is a directory, create the destination file name.
        if (destFile.type() == fileName::DIRECTORY)
        {
            destFile = destFile/src.name();
        }

        // Make sure the destination directory exists.
        if (!isDir(destFile.path()) && !mkDir(destFile.path()))
        {
            return false;
        }

        // Open and check streams.
        std::ifstream srcStream(src.c_str());
        if (!srcStream)
        {
            return false;
        }

        std::ofstream destStream(destFile.c_str());
        if (!destStream)
        {
            return false;
        }

        // Copy character data.
        char ch;
        while (srcStream.get(ch))
        {
            destStream.put(ch);
        }

        // Final check.
        if (!srcStream.eof() || !destStream)
        {
            return false;
        }
    }
    else if (src.type() == fileName::DIRECTORY)
    {
        // If dest is a directory, create the destination file name.
        if (destFile.type() == fileName::DIRECTORY)
        {
            destFile = destFile/src.component(src.components().size() -1);
        }

        // Make sure the destination directory exists.
        if (!isDir(destFile) && !mkDir(destFile))
        {
            return false;
        }

        // Copy files
        fileNameList contents = readDir(src, fileName::FILE, false);
        forAll(contents, i)
        {
            if (POSIX::debug)
            {
                Info<< "Copying : " << src/contents[i]
                    << " to " << destFile/contents[i] << endl;
            }

            // File to file.
            cp(src/contents[i], destFile/contents[i]);
        }

        // Copy sub directories.
        fileNameList subdirs = readDir(src, fileName::DIRECTORY);
        forAll(subdirs, i)
        {
            if (POSIX::debug)
            {
                Info<< "Copying : " << src/subdirs[i]
                    << " to " << destFile << endl;
            }

            // Dir to Dir.
            cp(src/subdirs[i], destFile);
        }
    }

    return true;
}


// Create a softlink. dst should not exist. Returns true if successful.
bool Foam::ln(const fileName& src, const fileName& dst)
{
    if (POSIX::debug)
    {
        Info<< "Create softlink from : " << src << " to " << dst
            << endl;
    }

    if (exists(dst))
    {
        WarningIn("ln(const fileName&, const fileName&)")
            << "destination " << dst << " already exists. Not linking."
            << endl;
        return false;
    }

    if (!exists(src))
    {
        WarningIn("ln(const fileName&, const fileName&)")
            << "source " << src << " does not exist." << endl;
        return false;
    }

    if (symlink(src.c_str(), dst.c_str()) == 0)
    {
        return true;
    }
    else
    {
        WarningIn("ln(const fileName&, const fileName&)")
            << "symlink from " << src << " to " << dst << " failed." << endl;
        return false;
    }
}


// Rename srcFile dstFile
bool Foam::mv(const fileName& src, const fileName& dst)
{
    if (POSIX::debug)
    {
        Info<< "Move : " << src << " to " << dst << endl;
    }

    if
    (
        dst.type() == fileName::DIRECTORY
     && src.type() != fileName::DIRECTORY
    )
    {
        const fileName dstName(dst/src.name());

        return rename(src.c_str(), dstName.c_str()) == 0;
    }
    else
    {
        return rename(src.c_str(), dst.c_str()) == 0;
    }
}


//- Rename to a corresponding backup file
//  If the backup file already exists, attempt with "01" .. "99" index
bool Foam::mvBak(const fileName& src, const std::string& ext)
{
    if (POSIX::debug)
    {
        Info<< "mvBak : " << src << " to extension " << ext << endl;
    }

    if (exists(src, false))
    {
        const int maxIndex = 99;
        char index[3];

        for (int n = 0; n <= maxIndex; n++)
        {
            fileName dstName(src + "." + ext);
            if (n)
            {
                sprintf(index, "%02d", n);
                dstName += index;
            }

            // avoid overwriting existing files, except for the last
            // possible index where we have no choice
            if (!exists(dstName, false) || n == maxIndex)
            {
                return rename(src.c_str(), dstName.c_str()) == 0;
            }

        }
    }

    // fall-through: nothing to do
    return false;
}



// Remove a file, returning true if successful otherwise false
bool Foam::rm(const fileName& file)
{
    if (POSIX::debug)
    {
        Info<< "Removing : " << file << endl;
    }

    // Try returning plain file name; if not there, try with .gz
    if (remove(file.c_str()) == 0)
    {
        return true;
    }
    else
    {
        return remove(string(file + ".gz").c_str()) == 0;
    }
}


// Remove a dirctory and its contents
bool Foam::rmDir(const fileName& directory)
{
    if (POSIX::debug)
    {
        Info<< "rmDir(const fileName&) : "
            << "removing directory " << directory << endl;
    }

    // Pointers to the directory entries
    DIR *source;
    struct dirent *list;

    // Attempt to open directory and set the structure pointer
    if ((source = opendir(directory.c_str())) == NULL)
    {
        WarningIn("rmDir(const fileName&)")
            << "cannot open directory " << directory << endl;

        return false;
    }
    else
    {
        // Read and parse all the entries in the directory
        while ((list = readdir(source)) != NULL)
        {
            fileName fName(list->d_name);

            if (fName != "." && fName != "..")
            {
                fileName path = directory/fName;

                if (path.type() == fileName::DIRECTORY)
                {
                    if (!rmDir(path))
                    {
                        WarningIn("rmDir(const fileName&)")
                            << "failed to remove directory " << fName
                            << " while removing directory " << directory
                            << endl;

                        closedir(source);

                        return false;
                    }
                }
                else
                {
                    if (!rm(path))
                    {
                        WarningIn("rmDir(const fileName&)")
                            << "failed to remove file " << fName
                            << " while removing directory " << directory
                            << endl;

                        closedir(source);

                        return false;
                    }
                }
            }

        }

        if (!rm(directory))
        {
            WarningIn("rmDir(const fileName&)")
                << "failed to remove directory " << directory << endl;

            closedir(source);

            return false;
        }

        closedir(source);

        return true;
    }
}


unsigned int Foam::sleep(const unsigned int s)
{
    return ::sleep(s);
}


void Foam::fdClose(const int fd)
{
    if (close(fd) != 0)
    {
        FatalErrorIn
        (
            "fdClose(const int fd)"
        )   << "close error on " << fd << endl
            << abort(FatalError);
    }
}


bool Foam::ping
(
    const word& destName,
    const label destPort,
    const label timeOut
)
{
    char *serverAddress;
    struct in_addr *ptr;
    struct hostent *hostPtr;
    volatile int sockfd;
    struct sockaddr_in destAddr;      // will hold the destination addr
    u_int addr;

    if ((hostPtr = gethostbyname(destName.c_str())) == NULL)
    {
        FatalErrorIn
        (
            "Foam::ping(const word&, const label)"
        )   << "gethostbyname error " << h_errno << " for host " << destName
            << abort(FatalError);
    }

    // Get first of the SLL of addresses
    serverAddress = *(hostPtr->h_addr_list);
    ptr = reinterpret_cast<struct in_addr*>(serverAddress);
    addr = ptr->s_addr;

    // Allocate socket
    sockfd = socket(AF_INET, SOCK_STREAM, 0);
    if (sockfd < 0)
    {
        FatalErrorIn
        (
            "Foam::ping(const word&, const label)"
        )   << "socket error"
            << abort(FatalError);
    }

    // Fill sockaddr_in structure with dest address and port
    memset (reinterpret_cast<char *>(&destAddr), '\0', sizeof(destAddr));
    destAddr.sin_family = AF_INET;
    destAddr.sin_port = htons(ushort(destPort));
    destAddr.sin_addr.s_addr = addr;


    timer myTimer(timeOut);

    if (timedOut(myTimer))
    {
        // Setjmp from timer jumps back to here
        fdClose(sockfd);
        return false;
    }

    if
    (
        connect
        (
            sockfd,
            reinterpret_cast<struct sockaddr*>(&destAddr),
            sizeof(struct sockaddr)
        ) != 0
    )
    {
        // Connection refused. Check if network was actually used or not.

        int connectErr = errno;

        fdClose(sockfd);

        if (connectErr == ECONNREFUSED)
        {
            return true;
        }
        //perror("connect");

        return false;
    }

    fdClose(sockfd);

    return true;
}


bool Foam::ping(const word& hostname, const label timeOut)
{
    return ping(hostname, 222, timeOut) || ping(hostname, 22, timeOut);
}


int Foam::system(const string& command)
{
    return ::system(command.c_str());
}


// ************************************************************************* //
