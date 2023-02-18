#!/usr/bin/env perl
use strict;
use warnings;

use AnyEvent;
use Linux::Inotify2;
use Getopt::Long qw(:config no_ignore_case);
use Pod::Usage;
use Config::Tiny;
use File::Find;
use Array::Utils;
use POSIX;
use File::Spec::Functions qw( catdir catfile rel2abs );
use Parallel::ForkManager;
use Mojo::SQLite;
use Mojo::Log;
use Fcntl qw/:flock/;

use Data::Printer;

our $VERSION = '0.1.0';

my %Watcher = ();
my %Seqpath = ();
my %Default = (
	pidfile => "/var/run/seqwatch.pid",
	logfile => "/var/log/seqwatch.log",
	confile => "/etc/seqwatch.conf",
	sqlite  => "",
	interval => 5,
	daemon  => 0,
	multiproc => 5,
	rundir => '.',
	fastqdir => '.',
	runrep => qr/\d{6,8}_[a-zA-Z\d\-]+_\d{1,}_[a-z-A-Z\d\-]+/,
	required => ['RunInfo.xml', 'RTAComplete.txt'],
	inchange  => 60,
	inprocess => 86400,
	bcl2fastq => 'bcl2fastq',
);
my @Need = qw/rundir fastqdir runrep required inchange inprocess/;

my %cmdline = ();
GetOptions(\%cmdline,
	'confile|c=s',
	'daemon|d',
	'sqlite|s=s',
	'pidfile|p=s',
	'logfile|l=s',
	'multiproc|m=i',
	'interval|t=i',
	'inchage|i',
	'inprocess|I=i',
	'required|r=s@',
	'runrep|e=s',
	'bcl2fastq|b',
	'policy|P=s',
	'help|h',
    'verbose|v',
) or pod2usage(2);
pod2usage(1) if delete $cmdline{help};

my $Verbose = delete $cmdline{verbose} || 0;
my $Policy = delete $cmdline{policy} || 'scan';

if(defined($cmdline{confile})){
	my $conf = Config::Tiny->new->read( $cmdline{confile} );
	my $default = delete $conf->{_} || {};
	%Default = (%Default, %$default);
	
	%Seqpath = %{$conf};
}

if(@ARGV){
	$Default{rundir} = shift;
	$Default{fastqdir} = shift || $Default{rundir};
}
	
%Default = (%Default, %cmdline);
$Seqpath{default} = \%Default;
foreach my $section (keys %Seqpath){
	die "'rundir' is required in section" unless exists $Seqpath{$section}{rundir};
	$Seqpath{$section}{fastqdir} ||= $Seqpath{$section}{rundir};
	$Seqpath{$section}{name} = $section;
	unless( -d $Seqpath{$section}{rundir}){
		die "dir no exists: $Seqpath{$section}{rundir}\n";
	}
	unless( -d $Seqpath{$section}{fastqdir} ){
		die "dir no exists: $Seqpath{$section}{fastqdir}\n";
	}
	foreach my $need ( @Need ){
		$Seqpath{$section}{$need} //= $Seqpath{default}{$need};
	}
    if($Seqpath{$section}{runrep} and ref $Seqpath{$section}{runrep} ne 'Regexp'){
        $Seqpath{$section}{runrep} = qr/$Seqpath{$section}{runrep}/;
    }
	#my $exceed_time = now_time(time() - $Seqpath{$section}{inchange});
	#$Seqpath{$section}{condsql} = {
	#	modify_time => { '<', $Seqpath{$section}{inchange} },
	#	begin_bcl2fastq => [ {'=', undef}, {'<', now_time(time() - $Seqpath{$section}{inprocess})} ],
	#	end_bcl2fastq => undef
	#};
}

#print Dumper(\%Seqpath);

my $Logger = Mojo::Log->new;
$Logger->handle( *STDERR );
if(open my $fh, ">$Default{logfile}"){
	$Logger->handle( $fh );
}

if($Default{daemon}){
	my $childpid = fork;

	if($childpid){
		open(FH, ">", "$Default{pidfile}") or pod2usage("can't open pidfile: $Default{pidfile}");
		print FH "$childpid";
        flock FH, LOCK_EX or die "daemon already exists!\n";
		#close(FH); 
	}

	die "Can't fork: $!\n" if(! defined $childpid);
  	exit if($childpid);

  	POSIX::setsid() or die "Can't start a new session: $!";
  	open STDIN, "</dev/null";
	open STDOUT, ">/dev/null";
	open STDERR, ">&STDOUT";
	umask 0;
    foreach my $section (keys %Seqpath){
    	$Seqpath{$section}{rundir} = rel2abs($Seqpath{$section}{rundir}) if( exists($Seqpath{$section}{rundir}) );
    }
#	chdir "/";
}

my $Sbatch = 0;
unless(system("sbatch --version")){
    $Sbatch = 1;
}
if(system("$Default{bcl2fastq} --version")){
	die "please -b bcl2fastq path\n";
}

my $Sqlite = '';
my $Inotify = new Linux::Inotify2;

if($Default{sqlite}){
	$Sqlite = Mojo::SQLite->new($Default{sqlite});
	$Sqlite->migrations->name('runs')->from_string(<<EOF)->migrate;
-- 1 up
CREATE TABLE runs (id INTEGER PRIMARY KEY AUTOINCREMENT, run_name TEXT NOT NULL UNIQUE, section TEXT NOT NULL, create_time DATETIME NOT NULL, modify_time DATETIME DEFAULT NULL, begin_bcl2fastq DATETIME DEFAULT NULL, end_bcl2fastq DATETIME DEFAULT NULL,files INT, dirs INT, state BOOLEAN);
-- 1 down
DROP TABLE names;
EOF
}

my ($i, $t);
foreach my $section (keys %Seqpath){
	if($Policy eq 'inotify'){
		pathDispatch(%{$Seqpath{$section}});
		$i = AnyEvent->io( fh => $Inotify->fileno, poll => 'r', cb => sub { $Inotify->poll } );
	}else{
		my %watch = %{$Seqpath{$section}};
		my $section = $watch{name};
	    foreach my $key (@Need){
	        my $info = ref $watch{$key} eq ref [] ? join ", ", @{$watch{$key}} : $watch{$key};
	        $Logger->info("[corak] scan $key: $info");
	    }

		$i = AnyEvent->timer( interval => $Default{interval}, cb => sub { pathScan( %watch ) });
	}
}

$t = AnyEvent->timer( interval => $Default{interval}, cb => \&checkCond);

AnyEvent->condvar->recv;

sub pathScan{
	my %watch = @_;

	my $section = $watch{name};
    my @current_dirs = map { s/$watch{rundir}\/*//; $_ } grep { -d and /$watch{runrep}/ } glob("$watch{rundir}/*");
    my @expire_dirs = map { s/$watch{fastqdir}\/*//; $_ } grep { -d and /$watch{runrep}/ and -f catfile($_, 'log') } glob("$watch{fastqdir}/*");
    my @current_runs = Array::Utils::array_minus(@current_dirs, @expire_dirs);
    my @old_runs = keys %{$Watcher{$section}};
    my @new_runs = Array::Utils::array_minus(@current_runs, @old_runs);

    # old run 
    for my $run (@old_runs){
    	my $path = catdir($watch{rundir}, $run);
    	unless($Watcher{$section}{$run}{begin_bcl2fastq}){
    		my ($files, $dirs) = (0, 0);
    		my $max_mtime = $Watcher{$section}{$run}{modify_time} || 0;
    		File::Find::find({wanted => sub {
    				$dirs++ if -d $File::Find::name;
    				if(-f $File::Find::name){
    					$files++;
    					my $mtime = [stat(_)]->[9];
    					$Logger->info("[change] " . $File::Find::name) if $Verbose and $Watcher{$section}{$run}{modify_time} and $mtime > $Watcher{$section}{$run}{modify_time};
    					$max_mtime = $mtime if $mtime > $max_mtime;
    				} 
    			}, no_chdir => 1,
    		}, $path);
    		$Watcher{$section}{$run}{files} = $files;
    		$Watcher{$section}{$run}{dirs} = $dirs;
    		$Watcher{$section}{$run}{modify_time} = $max_mtime;
    	}
    }
    # new run
    for my $run (@new_runs){
    	my $path = catdir($watch{rundir}, $run);
    	my $create_time = [stat($path)]->[9];
    	$Watcher{$section}{$run}{create_time} = $create_time;
    	$Logger->info("[create] $run");
    }
}

sub pathWatch{
	my %watch = @_;

	my $section = $watch{name};
    foreach my $key (@Need){
        my $info = ref $watch{$key} eq ref [] ? join ", ", @{$watch{$key}} : $watch{$key};
        $Logger->info("[corak] watch $key: $info");
    }
	$Inotify->watch($watch{rundir}, IN_CREATE, sub{
		my $e = shift;
		my $runpath = $e->fullname;
		my $mask = $e->mask;
		
		if($e->IN_ISDIR and $e->IN_CREATE and -d $runpath and $runpath =~ /$watch{runrep}/){
			$runpath =~ /($watch{runrep})/;
			my $run_name = $1;
			$Logger->info("[create] $run_name");
			$Watcher{$section}{$run_name}{create_time} = time;
			my $w = $Inotify->watch($runpath, IN_MODIFY|IN_CREATE, sub {
					my $e = shift;
					runChange($e, %watch);
				});
			unless($w){
				$Logger->falta("Inotify error: $!");
				die "can not create watcher, please see detail in log file\n";
			}
			push @{$Watcher{$section}{$run_name}{ws}}, $w;
		}
	});
}

sub runChange{
	my ($e, %watch) = @_;
	my $filename = $e->fullname;
	my $mask = $e->mask;

	$filename =~ /$watch{rundir}\/+($watch{runrep})\/*/;
	my $run_name = $1;
	my $section = $watch{name};
	$Logger->info("[change] $filename: " . (sprintf "%X", $mask)) if $Verbose;

	$Watcher{$section}{$run_name}{modify_time} = time;
	# if new dir, watch it 
	if($e->IN_CREATE){
        if(-d $filename){
		    my $w = $Inotify->watch($filename, IN_MODIFY|IN_CREATE, sub{
				my $e = shift;
				runChange($e, %watch);
			});
		    unless($w){
				$Logger->fatal("Inotify error: $!");
				die "can not create watcher, please see detail in log file\n";
			}
		    push @{$Watcher{$section}{$run_name}{ws}}, $w;
            $Watcher{$section}{$run_name}{dirs} += 1;         # 与实际数量少, 捕获丢失？
        }elsif(-f $filename){
            $Watcher{$section}{$run_name}{files} += 1;        # 与实际数量少，捕获丢失？
        }
	}

}

sub checkCond{
    #p %Watcher;
	foreach my $section (keys %Seqpath){
		my $inchange = $Seqpath{$section}{inchange};
		my $inprocess = $Seqpath{$section}{inprocess};
		foreach my $run_name (keys %{$Watcher{$section}}){
		    my $run_dir = my $indir = catdir($Seqpath{$section}{rundir}, $run_name);
		    my $outdir = catdir($Seqpath{$section}{fastqdir}, $run_name);
			my $ok = catfile($outdir, 'ok');
            my $run = $Watcher{$section}{$run_name};
            if($run->{state}){
                if($Sqlite){
                    $Sqlite->db->insert('runs', {
                                        run_name => $run_name, 
                                        section  => $section,
                                        create_time => now_time($run->{create_time}),
                                        modify_time => now_time($run->{modify_time}),
                                        begin_bcl2fastq => now_time($run->{begin_bcl2fastq}),
                                        end_bcl2fastq => now_time($run->{end_bcl2fastq}),
                                        files => $run->{files},
                                        dirs  => $run->{dirs},
                                        stat  => $run->{stat},
                    });
                }
                if($run->{state} == 0){ # 失败发邮件通知？
                   
                }
                my $slurm_out = "slurm-". $run->{eid} . '.out';
                if( -f $slurm_out ){
                    my $out = `cat $slurm_out`;
                    unlink($slurm_out);
                }
                delete $Watcher{$section}{$run_name};
                        
                next;
            }
            if($run->{begin_bcl2fastq}){  # 已经bcl2fastq
                if( ! -f $ok ){  
                    $Logger->info("[finish] $run_name");
                    $run->{end_bcl2fastq} = time;
                    $run->{state} = 1;
                }elsif(time() - $run->{begin_bcl2fastq} > $inprocess){   # bcl2fastq 超时， 视为失败, 或尝试查找相关进程
                    $Logger->info("[failed] $run_name");
                    $run->{state} = 0;
                }
                next;
            }
			if($run->{modify_time} and time() - $run->{modify_time} > $inchange){ # 一段时间文件无改变
				unless(grep { ! -e $_ } map {catfile $run_dir, $_ } @{$Seqpath{$section}{required}}){ # check require files
					# cancel watcher
					$_->cancel for @{$run->{ws}};
					delete $run->{ws};

				    $run->{begin_bcl2fastq} = time;
					mkdir($outdir) unless(-d $outdir);
					open my $fh, ">$ok";
					close($ok);
                    $Logger->info("[bcl2fastq] $run_name");
					my $id = RunFastq($indir, $outdir);
                    $run->{eid} = $id;
				}
			}
		}
	}
}

sub RunFastq{
	my ($indir, $outdir) = @_;

	my $shell = "$outdir/run.sh";
	open FH, ">$shell";
	print FH <<EOF;
#!/bin/bash

bcl2fastq -R $indir -o $outdir > $outdir/log 2>&1 && rm $outdir/ok
EOF
    if($Sbatch){
		my $sl = `sbatch -c 8 $shell`;
        $sl =~ /Submitted batch job (\d+)/;
        return $1;
	}else{
		my $pid = fork;
		if($pid){
			my $w = AnyEvent->child(pid => $pid, cb => sub {
				$Logger->info('[corak] bcl2fastq exit:' . $outdir);
				});
            return $pid;
		}else{
			system("bash $shell");
		}
	}
}

sub now_time{
	my $time = shift || time;
	POSIX::strftime "%e-%b-%Y %H:%M:%S", localtime($time);
}

__END__

=encoding utf-8

=head1 NAME

seqwatch - 监控目录，如果下机数据传输完成，执行bcl2fastq

=head1 SYNOPSIS

seqwatch [options] run_dir fastq_dir

    run_dir:  下机数据上传目录

    fastq_dir: fastq文件生成目录
  
    Options:
         -h, --help                  命令帮助
         -P, --policy                scan or inotify
         -d, --daemon                以后台服务的方式运行
         -p, --pidfile               储存 PID 的文件
         -c, --confile               配置文件
         -s, --sqlite                Sqlite数据库文件
         -l, --logfile               日志文件
         -b, --bcl2fastq             bcl2fastq 命令路径
         -t, --interval              检测是否上传完毕的时间间隔
         -r, --required              Run中必须存在的文件
         -e，--runrep                Run目录名的正则表达式
         -i, --inchange              视为上传完毕的时间间隔
         -I, --inprocess             视为bcl2fastq完毕的时间间隔
         -v, --verbose               

    详情请执行：perldoc seqwatch

=head1 DESCRIPTION

监控二代测序下机数据是否已经完成上传到指定目录，如果完成，执行bcl2fastq生成fastq文件。监控有两种机制：一是
使用Linux的inotify功能，每当指定目录下产生一个创建或修改事件，inotify都会将事件传递给进程，进程对事件做相应
反应；它在分布式文件系统上（如PGPFS）会失效。 二是每隔一段时间扫描目录，以断定文件是否发生改变。

配置文件示例：

# example.conf

  logfile  = ./example.log
  sqlite   = ./test.db
  interval = 5
  runrep   = /\d{6,8}_[a-zA-Z\d\-]+_\d{1,}_[a-z-A-Z\d\-]+/

  [nextseq550]

  rundir = /share/Data/Test/nextseq550
  fastqdir = /share/Data/Test/nextseq550
  inchange = 60


  [hiseq2500]

  rundir = /share/Data/Test/nextseq550
  fastqdir = ./
  runrep = /^HISEQ_[a-zA-Z\d\-]+_\d{1,}_[a-z-A-Z\d\-]+/

=head1 AUTHOR

chenshulin@bdlifescience.com

=head1 BUGS

存在的问题：
   
   1. 因为inotify是借用内核的notify通知链技术，所以只有当运行该程序的机器内核对磁盘进行读写操作才会被监测，
   如果是另外一台机器对监控的目录进行读写，则不能监控到。

   2. 文件和目录统计不准确，原因模糊

   3. 没有一个充分的机制判断数据上传完成，仅仅根据断传的时间间隔太过武断; 将来可结合监控上传进程


