
# Initialize the environment with path variables, CFLAGS, and so on
bioinfo_path = '#lib/bioinfo-libs'
commons_path = '#lib/common-libs'
#math_path = '#libs/math'

vars = Variables('buildvars.py')

compiler = ARGUMENTS.get('compiler', 'gcc')

env = Environment(tools = ['default', 'packaging'],
      		  CC = compiler,
                  variables = vars,
                  CFLAGS = '-std=c99 -D_XOPEN_SOURCE=600 -D_GNU_SOURCE -fopenmp',
                  CPPPATH = ['#', '#src', '#include', bioinfo_path, commons_path, "%s/commons/argtable" % commons_path, "%s/commons/config" % commons_path ],
                  LIBPATH = [commons_path, bioinfo_path],
                  LIBS = ['m', 'z', 'curl'],
                  LINKFLAGS = ['-fopenmp'])

if int(ARGUMENTS.get('debug', '0')) == 1:
    debug = 1
    env['CFLAGS'] += ' -O0 -g'
else:
    debug = 0
    env['CFLAGS'] += ' -O3 -g'

env['objects'] = []

# Targets

SConscript(['%s/SConscript' % bioinfo_path,
            '%s/SConscript' % commons_path
            ], exports = ['env', 'debug', 'compiler'])

env.Program('hpg-bam',
             source = [Glob('src/*.c'),
                       "%s/libbioinfo.a" % bioinfo_path,
                       "%s/libcommon.a" % commons_path
                      ]
           )

env.Install('#bin', 'hpg-bam')

'''
if 'debian' in COMMAND_LINE_TARGETS:
    SConscript("deb/SConscript", exports = ['env'] )
'''
