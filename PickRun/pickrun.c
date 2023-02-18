/*
 * pickrun - run bcl2fastq when a run is uploaded
 */

/* 
 * Written by chenhsulin@bdlifescience.com
 */

#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>
#include <stdarg.h>
#include <errno.h>
#include <getopt.h>
#include <unistd.h>
#include <fcntl.h>
#include <signal.h>
#include <dirent.h>
#include <time.h>
#include <syslog.h>
#include <regex.h>
#include <ftw.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/resource.h>
/* max section watched */
#define MAX_SEC 3
#define MAX_REQ 4
#define MAX_PATH 1024
#define MAX_WORD 128
#define MAX_LINE 1024

static char const *program_name = "pickrun";
static char const *Version = "0.1.0";

//extern int log_to_stderr;
extern int errno;
extern char *optarg;
/*
 * Section is a set of run have same features
 */

typedef struct Run RUN;
typedef struct Section SECTION;

/* Section descriptor structure */

struct Section {
    char s_name[MAX_WORD];
    char s_rundir[MAX_PATH];
    char s_fqdir[MAX_PATH];
    unsigned short s_inchange;
    unsigned int s_inprocess;
    char *s_regex;
    unsigned s_num;
    RUN *s_head;
    char *s_required[MAX_REQ];
};


/* Run descriptor structure. */

enum STATE { CRT, UPD, PCS, FIS };

struct Run {
    char r_name[MAX_WORD];                     /* the run name, also the run's directory name */
    time_t r_mtime;                         /* last modify time iin run */
    enum STATE r_state;                     /* run state */
    SECTION *r_section;                      /* secton the run belong to */
    RUN *r_next;                            /* next run in section */
};

/* forward declarations */

static void msg_doit(int, int, const char *, va_list ap);
static void err_msg(const char *, ...);
static void log_msg(const char *, ...);
static void emit_try_help(void);
static void usage(int);
static void parse_conf(SECTION *, const char *);
static void pickrun(SECTION *);
static void debug(void);

static RUN* create_run(char *, time_t, enum STATE, SECTION *, RUN *); 
static void update_run(RUN *);
static void submit_run(RUN *);

static char* concat(const char *, ...);
static int files_exists(const char *, ...);
static int latest_time(const char *, const struct stat *, int);
static void daemonize(const char *cmd);

/* Option values. */

/* if true, daemonize (default false). */
static bool is_daemon;

/* if ture, use sqlite to store run info (default false). */
static bool is_sqlite;

/* if true, multi thread (default false). */
static bool is_thread;

/* if true, send mail when error (default false). */
static bool is_mail;

/* if true, use atime to determine run uploaded other than mtime (default false). */
static bool is_atime = 0;

/* default log to stderr */

static int log_to_stderr = 1;

/* iterval to scan new run (default INTERVAL). */
static short interval = 5;

/* pid file. */
static char *pidfile = "test.pid";

/* latest time */
static time_t latest;

/* default section */
static SECTION sections[MAX_SEC] = {
    {"default", ".", ".", 60, 86400, "\\d{6,8}_[a-zA-Z\\d\\-]+_\\d{1,}_[a-zA-Z\\d\\-]+", 0, NULL, {"RunInfo.xml", "SampleSheet.csv", "RTAComplete.txt"}}
};

/* options */
static struct option const longopts[] = {
    {"daemon", optional_argument, NULL, 'd'},
    {"sqlite", optional_argument, NULL, 's'},
    {"thread", optional_argument, NULL, 't'},
    {"mail",   optional_argument, NULL, 'm'},
    {"atime",  no_argument, NULL, 'a'},
    {"interval", required_argument, NULL, 'i'},
    {"pidfile", required_argument, NULL, 'p'},
    {"inchange", required_argument, NULL, 'I'}, 
    {"inprocess", required_argument, NULL, 'L'},
    {"regex", required_argument, NULL, 'r'},
    {"required", required_argument, NULL, 'R'},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'v'},
    {NULL, 0, NULL, 0}
};


int main(int argc, char *argv[])
{
    int optc;
    char *conf;

    while(( optc = getopt_long(argc, argv, "hd::s::t::m::ai:p:I:L:r:R;f:l:", longopts, NULL)) != -1 ){
        switch(optc){
        case 'd':
            is_daemon = 1;
            break;
        case 's':
            is_sqlite = 1;
            break;
        case 't':
            is_thread = 1;
            break;
        case 'm':
            is_mail = 1;
            break;
        case 'a':
            is_atime = 1;
        case 'f':
            conf = optarg;
            parse_conf(sections, conf);
            break;
        case 'h':
        case '?':
            usage(1);
            break;
        }
    }
    while(1)
    {
        pickrun(&sections[0]);
        sleep(5);
        debug();
    }
}


static void msg_doit(int errnoflag, int error, const char *fmt, va_list ap)
{
    char buf[MAX_LINE];

    vsnprintf(buf, MAX_LINE - 1, fmt, ap);
    if(errnoflag)
        snprintf(buf + strlen(buf), MAX_LINE - strlen(buf) - 1, ": %s", strerror(error));
    strcat(buf, "\n");
    if(log_to_stderr){
        fflush(stdout);
        fputs(buf, stderr);
        fflush(NULL);
    }else{
        syslog(1, "%s", buf);
   }
}

static void err_msg(const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    msg_doit(0, 0, fmt, ap);
    va_end(ap);
    exit(1);
}

static void log_msg(const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    msg_doit(0, 0, fmt, ap);
    va_end(ap);
}

static void emit_try_help(void)
{
    fprintf(stderr, "Try '%s --help' for more information.\n", program_name);
}

static void usage(int status)
{
    if(status == 0)
        emit_try_help();
    else{
        printf("\n\
pickrun version: %s\n\
USAGE: %s [OPTION] rundir fadir\n\
",
            Version, program_name);
        puts("\
Optional paremeters:\n\
    -d, --daemon [test.pid]         run as daemon\n\
    -s, --sqlite [test.db]          save info to sqlite\n\
        ");
    }
    exit(0);
}

static void parse_conf(SECTION *sec, const char *conf)
{
    int n, i, j, sec_idx;
    FILE *fh;
    char *tmp;
    char key[MAX_WORD], val[MAX_WORD];
    char line[MAX_LINE]; // pattern[MAX_LINE];

    if((fh = fopen(conf, "r")) == NULL)
        err_msg("fopen error: %s", conf);
    //sprintf(line, "%%%ds%%%ds", MAX_WORD-1, MAX_WORD-1);
    sec_idx = 0;
    while(fgets(line, MAX_LINE, fh) != NULL)
    {
        if(line[0] == '#' || isspace(line[0])) continue;
        if(line[0] == '[')
        {
            if(++sec_idx == MAX_SEC) break;
            memcpy(sec+sec_idx, sec, sizeof(struct Section));
            memset(sec[sec_idx].s_name, 0, MAX_WORD);
            for(i=1; line[i] != ']' && line[i] != '\0'; i++) sec[sec_idx].s_name[i-1] = line[i];
            continue;
        }
        memset(key, 0, sizeof(key));
        memset(val, 0, sizeof(val));
        n = sscanf(line, "%[^:]:%*[ \t]%s", key, val);
        //printf("line: %skey: %s\nval: %s\n", line, key, val);
        if(n == 2)
        {
            if(strcmp(key, "rundir") == 0)
                memcpy(sec[sec_idx].s_rundir, val, strlen(val) + 1);
            else if(strcmp(key, "fqdir") == 0)
                memcpy(sec[sec_idx].s_fqdir, val, strlen(val) + 1);
            else if(strcmp(key, "inchange") == 0)
                sec[sec_idx].s_inchange = atoi(val);
            else if(strcmp(key, "inprocess") == 0)
                sec[sec_idx].s_inprocess = atoi(val);
            else if(strcmp(key, "regex") == 0)
            {
                n = strlen(val);
                sec[sec_idx].s_regex = (char *)malloc(n+1);

                memcpy(sec[sec_idx].s_regex, val, n+1);
            }else if(strcmp(key, "required") == 0)
            {
                i = 0; j = 0; n = 0;
                tmp = (char *)malloc(MAX_WORD + 1);
                memset(tmp, 0, sizeof(tmp));
                while(val[i] != '\0')
                {
                    if(val[i] == ',')
                    {
                        sec[sec_idx].s_required[n++] = tmp;
                        tmp = (char *)malloc(MAX_WORD + 1);
                        memset(tmp, 0, sizeof(tmp));
                        j = 0; i++;
                        continue;
                    }
                    if(val[i] == ' ')
                        continue;
                    tmp[j] = val[i];
                    i++; j++;
                }
                sec[sec_idx].s_required[n] = tmp;
                while(n++ < MAX_SEC - 1)
                    sec[sec_idx].s_required[n] = NULL;
                    
            }
            else
                log_msg("Unkown key: %s", key);
        }
    }

}

static void pickrun(SECTION *sec)
{
    DIR *dp;
    struct dirent *dirp;
    struct stat statbuf;

    RUN *point_run, *tmp_run;
    
    if((dp = opendir(sec->s_rundir)) == NULL)
        err_msg("can't open directory: %s", sec->s_rundir);

    point_run = sec->s_head;
    while((dirp = readdir(dp)) != NULL){
        // filter directory . and ..
        if(strcmp(dirp->d_name, ".") == 0 || strcmp(dirp->d_name, "..") == 0)
            continue;
        
        memset(&statbuf, 0, sizeof(struct stat));
        if(lstat(dirp->d_name, &statbuf) == -1)
            err_msg(strerror(errno));
        
        // Isn't directory
        if((statbuf.st_mode & S_IFMT) != S_IFDIR)
            continue;
        // dir name don't match regex
        if(0){}
        
        /* it's a run, begin */
        // no run in list
        if(sec->s_head == NULL)
        {
            tmp_run = create_run("nil", 0, CRT, sec, NULL);
            sec->s_head = create_run(dirp->d_name, time(NULL), CRT, sec, tmp_run);
            point_run = sec->s_head;
        }
        //run already exists
        else if(strcmp(point_run->r_name, dirp->d_name) == 0)
        {
            switch(point_run->r_state)
            {
                case CRT:
                    update_run(point_run);
                    break;
                case UPD:
                    submit_run(point_run);
                    break;
                case PCS:
                    break;
                case FIS:
                    break;
                default:
                    break;
             }
        }
        // new run, insert to run list
        else
        {
            if((tmp_run = (RUN *)malloc(sizeof(RUN))) == NULL)
                err_msg(strerror(errno));
            memcpy(tmp_run, point_run, sizeof(RUN));
            memcpy(point_run->r_name,dirp->d_name, strlen(dirp->d_name) + 1);
            point_run->r_mtime = time(NULL);
            point_run->r_state = CRT;
            point_run->r_next = tmp_run;
        }
        point_run = point_run->r_next;
    }
}

static RUN *create_run(char *name, time_t mtime, enum STATE state, SECTION *sec, RUN *next)
{
    RUN *r;
    
    if((r = (RUN *)malloc(sizeof(RUN))) == NULL)
        err_msg(strerror(errno));

    memset(r, 0, sizeof(RUN));
    memcpy(r->r_name, name, strlen(name) + 1);
    r->r_mtime = mtime;
    r->r_state = state;
    r->r_section = sec;
    r->r_next = next;

    return(r);
}

static void update_run(RUN *r)
{
    time_t now;
    char *runbase;

    now = time(NULL);
    runbase = concat(r->r_section->s_rundir, "/", r->r_name, NULL);

    // 下机项目中的文件长时间没有变化，可能下机完成
    if(now - r->r_mtime > r->r_section->s_inchange)
    {
        int all_file_exists;
        char **filename;
        char *path;
        
        all_file_exists = 1;
        for(filename = r->r_section->s_required; *filename != NULL; filename++)
        {
            path = concat(runbase, "/", *filename, NULL);
            printf("required: %s\n", path);
            if(access(path,F_OK) == -1)
                all_file_exists = 0;
            free(path);
        }
        if(all_file_exists)
            r->r_state = UPD;
    }
    // 正在上传下机结果
    else
    {
        latest = r->r_mtime;
        if(ftw(runbase, latest_time, 10) == -1)
            log_msg("scan run directory: %s", strerror(errno));
        r->r_mtime = r->r_mtime > latest ? r->r_mtime : latest;
    }
}

static void submit_run(RUN *r)
{
    char *indir, *outdir;
    char *argv[10];
    int i, fd;
    pid_t pid;

    i = 0;
    indir = concat(r->r_section->s_rundir, "/", r->r_name, '\0');
    outdir = concat(r->r_section->s_fqdir, "/", r->r_name, '\0');
    if(execlp("srun", "srun", (char *)0) < 0)  // no slurm
    {
        argv[i++] = "bcl2fastq";
    }
    else
    {
        argv[i++] = "srun";
        argv[i++] = "bcl2fastq";
    }

    argv[i++] = "--no-lane-splitting";

    argv[i++] = "-R";
    argv[i++] = indir;

    argv[i++] = "-o";
    argv[i++] = outdir;

    argv[i++] = "--interop-dir";
    argv[i++] = outdir;

    while(i<10)
        argv[i++] = NULL;

    char *check_file = "ok";
    if((fd = create(check_file, FILE_MODE)) < 0)
        err_msg("can't create file: %s", check_file);
    else
    {
        if((pid = vfork()) < 0)
            err_log("vfork error!")
        else if(pid == 0)     /* child */
        {
            char *log_file = "log";
            int fd2;


            if(dup2(fd2, 1) < 0)
                err_msg("dup2 for stdout error");
            if(dup2(fd2, 2) < 0)
                err_msg("dup2 for stferr error");

            execvp(argv[0], argv);
            _exit(127);
        }
        else
        {
            close(fd);
            return 0;
        }
   }
}
    

static void debug(void)
{
    int i, j;

    for(i=0; i<MAX_SEC; i++)
    {
        if(strcmp(sections[i].s_name, "") != 0)
        {
            printf("[%s]\nrundir: %s\nfqdir: %s\ninchange: %d\ninprocess: %d\nregex: %s\n",sections[i].s_name,
                    sections[i].s_rundir, sections[i].s_fqdir, sections[i].s_inchange, sections[i].s_inprocess,
                    sections[i].s_regex);
            
	        printf("required: ");
                for(j=0; j<MAX_REQ; j++)
                    if(sections[i].s_required[j] != NULL && strcmp(sections[i].s_required[j], "") != 0)   // C 中字符串数组未初始化的元素被初始化为NULL，而不是“”
                        printf("%s ", sections[i].s_required[j]);
            
            if(sections[i].s_head != NULL)
            {
                RUN *r;
                
                printf("\n");
                for(r = sections[i].s_head; r != NULL; r = r->r_next)
                    printf("%s:%ld:%ld\n", r->r_name, r->r_mtime, r->r_state);
                printf("##################################################################\n\n");
            }
            else
                printf("\n");
        }
    }

}

static char* concat(const char *str, ...)
{
    va_list ap;
    size_t allocated = 100;
    char *result = (char *)malloc(allocated);

    if(result != NULL)
    {
        char *newp;
        char *wp;
        const char *s;

        va_start(ap, str);

        wp = result;
        for(s = str; s != NULL; s = va_arg(ap, const char *))
        {
            size_t len = strlen(s);

            if(wp + len + 1 > result + allocated)
            {
                allocated = (allocated + len) * 2;
                newp = (char *)realloc(result, allocated);
                if(newp == NULL)
                {
                    free(result);
                    return NULL;
                }
                wp = newp + (wp - result);
                result = newp;
            }

            wp = mempcpy(wp, s, len);
        }

        *wp++ = '\0';

        newp = realloc(result, wp - result);
        if(newp != NULL)
            result = newp;
        va_end(ap);
    }
    return result;
}

static int files_exists(const char *path, ...)
{
    va_list ap;
    int all_exists = 1;
    const char *p;

    va_start(ap, path);

    for(p = path; p != NULL; p = va_arg(ap, const char *))
        if(access(p, F_OK) == -1)
            all_exists = 0;

    return all_exists;
}

static int latest_time(const char *fpath, const struct stat *statbuf, int typeflag)
{
    if(typeflag == FTW_F)
        if(is_atime)
        {
            if(statbuf->st_atime > latest)
                latest = statbuf->st_atime;
        }
        else
        {
            if(statbuf->st_mtime > latest)
                latest = statbuf->st_mtime;
        }
    //printf("%s:%ld:%ld\n", fpath, latest, statbuf->st_mtime);
    return 0;
}

static void daemonize(const char *cmd)
{
    int i, fd0, fd1, fd2;
    pid_t pid;
    struct rlimit rl;
    struct sigaction sa;

    /*
     * Clear file creation mask
     */
    umask(0);
    if(getrlimit(RLIMIT_NOFILE, &rl) < 0)
        err_msg("%s: can't get file limit", cmd);

    /*
     * Become a session leader to lose controlling TTY.
     */
    if((pid = fork()) < 0)
        err_msg("%s: can't fork", cmd);
    else if(pid != 0)
        exit(0);
    setsid();

    /*
     * Ensuer future opens won't allocate controlling TTYs.
     */
    sa.sa_handler = SIG_IGN;
    sigemptyset(&sa.sa_mask);
    sa.sa_flags = 0;
    if(sigaction(SIGHUP, &sa, NULL) < 0)
    {
        fprintf(stderr, "%s: can't ignore SIGHUP", cmd);
        exit(255);
    }
    if((pid = fork()) < 0)
    {
        fprintf(stderr, "%s: can't fork", cmd);
        exit(255);
    }
    else if(pid != 0)
        exit(0);

    /*
     * Change the current working directory to the root so 
     * we won't prevent filesystems from being unmount.
     */
    if(chdir("/") < 0)
    {
        fprintf(stderr, "%s: can't change directory to /", cmd);
        exit(255);
    }

    /*
     * Close all open file descriptors.
     */
    if(rl.rlim_max == RLIM_INFINITY)
        rl.rlim_max = 1024;
    for(i=0; i<rl.rlim_max; i++)
        close(i);

    /* 
     * Attach file descriptors 0, 1, and 2 to /dev/null.
     */
    fd0 = open("/dev/null", O_RDWR);
    fd1 = dup(0);
    fd2 = dup(0);

    /*
     * Initialize the log file.
     */
    openlog(cmd, LOG_CONS, LOG_DAEMON);
    if(fd0 != 0 || fd1 != 1 || fd2 !=2)
    {
        syslog(LOG_ERR, "unexpected file descriptors %d %d %d", fd0, fd1, fd2);
        exit(1);
    }
}
