#! /usr/bin/env python
import sys, os
import os.path
import tempfile
import urllib
import shutil
import subprocess
import atexit

class OfflineConverter:

    databases = {}
    databases['orcoff'] = ( '-t', 'oracle', '-h', 'cmsr1-s.cern.ch', '-d', 'cms_cond.cern.ch', '-u', 'cms_hlt_gui_r',     '-s', 'convertme!' )
    databases['hltdev'] = ( '-t', 'oracle', '-h', 'cmsr1-s.cern.ch', '-d', 'cms_cond.cern.ch', '-u', 'cms_hltdev_reader', '-s', 'convertme!' )


    @staticmethod
    def CheckTempDirectory(dir):
        dir = os.path.realpath(dir)
        if not os.path.isdir(dir):
            try:
                os.makedirs(dir)
            except:
                return None
        return dir


    def __init__(self, database = 'hltdev', url = None, verbose = False):
        self.verbose = verbose
        self.baseDir = '/afs/cern.ch/user/c/confdb/www/lib'
        self.baseUrl = 'http://confdb.web.cern.ch/confdb/lib'
        self.jars = ( 'ojdbc6.jar', 'cmssw-evf-confdb-converter.jar' )
        self.workDir = ''

        # check the database
        if database in self.databases:
            # load the connection parameters for the given database
            self.connect = self.databases[database]
        else:
            # unsupported database
            sys.stderr.write( "ERROR: unknown database \"%s\"\n" % database)
            sys.exit(1)

        # check for a custom base URL
        if url is not None:
            self.baseUrl = url

        # try to read the .jar files from AFS, or download them
        if os.path.isdir(self.baseDir) and all(os.path.isfile(self.baseDir + '/' + jar) for jar in self.jars):
            # read the .jar fles from AFS
            self.workDir = self.baseDir
        else:
            # try to use $CMSSW_BASE/tmp
            self.workDir = OfflineConverter.CheckTempDirectory(os.environ['CMSSW_BASE'] + '/tmp/confdb')
            if not self.workDir:
                # try to use $TMP
                self.workDir = OfflineConverter.CheckTempDirectory(os.environ['TMP'] + '/confdb')
            if not self.workDir:
                # create a new temporary directory, and install a cleanup callback
                self.workDir = tempfile.mkdtemp()
                atexit.register(shutil.rmtree, self.workDir)
            # download the .jar files
            for jar in self.jars:
                # check if the file is already present
                if os.path.exists(self.workDir + '/' + jar):
                    continue
                # download to a temporay name and use an atomic rename (in case an other istance is downloading the same file
                handle, temp = tempfile.mkstemp(dir = self.workDir, prefix = jar + '.')
                os.close(handle)
                urllib.urlretrieve(self.baseUrl + '/' + jar, temp)
                if not os.path.exists(self.workDir + '/' + jar):
                    os.rename(temp, self.workDir + '/' + jar)
                else:
                    os.unlink(temp)

        # setup the java command line and CLASSPATH
        if self.verbose:
            sys.stderr.write("workDir = %s\n" % self.workDir)
# Use non-blocking random # source /dev/urandom (instead of /dev/random), see:
# http://blockdump.blogspot.fr/2012/07/connection-problems-inbound-connection.html
# Also deal with timezone region not found
# http://stackoverflow.com/questions/9156379/ora-01882-timezone-region-not-found
        self.javaCmd = ( 'java', '-cp', ':'.join(self.workDir + '/' + jar for jar in self.jars),'-Djava.security.egd=file:///dev/urandom','-Doracle.jdbc.timezoneAsRegion=false','confdb.converter.BrowserConverter' )


    def query(self, *args):
        args = self.javaCmd + self.connect + args
        if self.verbose:
            sys.stderr.write("\n" + ' '.join(args) + "\n\n" )
        sub = subprocess.Popen(
            args,
            stdin = None,
            stdout = subprocess.PIPE,
            stderr = subprocess.PIPE,
            shell = False,
            universal_newlines = True )
        return sub.communicate()

def help():
    sys.stdout.write("""Usage: %s OPTIONS

        --hltdev|--orcoff           (target db [default: hltdev])

        --configId <id>             (specify configuration by id)
        --configName <name>         (specify configuration by name)
        --runNumber <run>           (specify configuration by run)
          [exactly one of --configId OR --configName OR --runNumber is required]

        --cff                       (retrieve configuration *fragment*)
        --input <f1.root[,f2.root]> (insert PoolSource with specified fileNames)
        --input <files.list>        (read a text file which lists input ROOT files)
        --output <out.root>         (insert PoolOutputModule w/ specified fileName)
        --nopsets                   (exclude all globale psets)
        --noedsources               (exclude all edsources)
        --noes                      (exclude all essources *and* esmodules)
        --noessources               (exclude all essources)
        --noesmodules               (exclude all esmodules)
        --noservices                (exclude all services)
        --nooutput                  (exclude all output modules)
        --nopaths                   (exclude all paths [+=referenced seqs&mods])
        --nosequences               (don't define sequences [+=referenced s&m])
        --nomodules                 (don't define modules)
        --psets <pset1[,pset2]>     (include only specified global psets)
        --psets <-pset1[,-pset2]>   (include all global psets but the specified)
        --essources <ess1[,ess2]>   (include only specified essources)
        --essources <-ess1[,-ess2]> (include all essources but the specified)
        --esmodules <esm1[,esm2]>   (include only specified esmodules)
        --esmodules <-esm1[,-esm2]> (include all esmodules but the specified)
        --services <svc1[,svc2]>    (include only specified services)
        --services <-svc1[,-svc2]>  (include all services but the specified)
        --paths <p1[,p2]>           (include only specified paths)
        --paths <-p1[,-p2]>         (include all paths but the specified)
        --streams <s1[,s2]>         (include only specified streams)
        --datasets <d1[,d2]>        (include only specified datasets)
        --sequences <s1[,s2]>       (include sequences, referenced or not!)
        --modules <p1[,p2]>         (include modules, referenced or not!)
        --blocks <m1::p1[,p2][,m2]> (generate parameter blocks)
""")


def main():
    args = sys.argv[1:]
    db = 'hltdev'

    if not args:
        help()
        sys.exit(1)

    if '--help' in args or '-h' in args:
        help()
        sys.exit(0)

    if '--orcoff' in args and '--hltdev' in args:
        sys.stderr.write( "ERROR: conflicting database specifications \"--hltdev\" and \"--orcoff\"\n" )
        sys.exit(1)

    if '--runNumber' in args and '--hltdev' in args:
        sys.stderr.write( "ERROR: conflicting database specifications \"--hltdev\" and \"--runNumber\"\n" )
        sys.exit(1)

    if '--hltdev' in args:
        db = 'hltdev'
        args.remove('--hltdev')

    if '--orcoff' in args:
        db = 'orcoff'
        args.remove('--orcoff')

    if '--runNumber' in args:
        db = 'orcoff'

    converter = OfflineConverter(database = db)
    out, err = converter.query( * args )
    if 'ERROR' in err:
        sys.stderr.write( "%s: error while retriving the HLT menu\n\n%s\n\n" % (sys.argv[0], err) )
        sys.exit(1)
    else:
        sys.stdout.write( out )


if __name__ == "__main__":
    main()
