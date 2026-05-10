#include <R.h>
#include <Rinternals.h>

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include <pandaseq.h>
#include <pandaseq-args.h>

typedef struct {
    int saved_stderr_fd;
    FILE *fp;
} stderr_capture_t;

static char *xstrdup_local(const char *s)
{
    size_t n;
    char *out;

    if (s == NULL) {
        return NULL;
    }

    n = strlen(s);
    out = (char *) calloc(n + 1, 1);
    if (out == NULL) {
        return NULL;
    }

    memcpy(out, s, n);
    out[n] = '\0';
    return out;
}

static char *read_all_from_fp(FILE *fp)
{
    long size;
    char *buf;
    size_t nread;

    if (fp == NULL) {
        return xstrdup_local("");
    }

    if (fflush(fp) != 0) {
        return xstrdup_local("");
    }

    if (fseek(fp, 0, SEEK_END) != 0) {
        return xstrdup_local("");
    }

    size = ftell(fp);
    if (size < 0) {
        return xstrdup_local("");
    }

    if (fseek(fp, 0, SEEK_SET) != 0) {
        return xstrdup_local("");
    }

    buf = (char *) calloc((size_t) size + 1, 1);
    if (buf == NULL) {
        return xstrdup_local("");
    }

    nread = fread(buf, 1, (size_t) size, fp);
    buf[nread] = '\0';
    return buf;
}

static int stderr_capture_begin(stderr_capture_t *cap)
{
    if (cap == NULL) {
        return -1;
    }

    memset(cap, 0, sizeof(*cap));
    cap->saved_stderr_fd = -1;

    cap->fp = tmpfile();
    if (cap->fp == NULL) {
        return -1;
    }

    fflush(stderr);
    cap->saved_stderr_fd = dup(fileno(stderr));
    if (cap->saved_stderr_fd < 0) {
        fclose(cap->fp);
        cap->fp = NULL;
        return -1;
    }

    if (dup2(fileno(cap->fp), fileno(stderr)) < 0) {
        close(cap->saved_stderr_fd);
        cap->saved_stderr_fd = -1;
        fclose(cap->fp);
        cap->fp = NULL;
        return -1;
    }

    return 0;
}

static char *stderr_capture_end(stderr_capture_t *cap)
{
    char *buf;

    if (cap == NULL || cap->fp == NULL) {
        return xstrdup_local("");
    }

    fflush(stderr);

    if (cap->saved_stderr_fd >= 0) {
        dup2(cap->saved_stderr_fd, fileno(stderr));
        close(cap->saved_stderr_fd);
        cap->saved_stderr_fd = -1;
    }

    buf = read_all_from_fp(cap->fp);
    fclose(cap->fp);
    cap->fp = NULL;

    return buf;
}

static const char *find_option_value(int argc, char **argv, const char *opt)
{
    int i;

    for (i = 0; i < argc - 1; i++) {
        if (strcmp(argv[i], opt) == 0) {
            return argv[i + 1];
        }
    }

    return NULL;
}

static char *join_messages(const char *a, const char *b)
{
    size_t na = (a != NULL) ? strlen(a) : 0;
    size_t nb = (b != NULL) ? strlen(b) : 0;
    char *out;

    if (na == 0 && nb == 0) {
        return xstrdup_local("");
    }

    out = (char *) calloc(na + nb + 2, 1);
    if (out == NULL) {
        return xstrdup_local("");
    }

    if (na > 0) {
        memcpy(out, a, na);
    }

    if (na > 0 && nb > 0) {
        out[na] = '\n';
        memcpy(out + na + 1, b, nb);
    } else if (nb > 0) {
        memcpy(out, b, nb);
    }

    out[na + (na > 0 && nb > 0 ? 1 : 0) + nb] = '\0';
    return out;
}

static SEXP make_result(
    int ok,
    int status,
    const char *phase,
    const char *output_fastq,
    const char *log_file,
    const char *stderr_text
) {
    SEXP res, names;
    int i = 0;

    PROTECT(res = allocVector(VECSXP, 6));
    PROTECT(names = allocVector(STRSXP, 6));

    SET_STRING_ELT(names, 0, mkChar("ok"));
    SET_STRING_ELT(names, 1, mkChar("status"));
    SET_STRING_ELT(names, 2, mkChar("phase"));
    SET_STRING_ELT(names, 3, mkChar("output_fastq"));
    SET_STRING_ELT(names, 4, mkChar("log_file"));
    SET_STRING_ELT(names, 5, mkChar("stderr"));

    SET_VECTOR_ELT(res, 0, ScalarLogical(ok));
    SET_VECTOR_ELT(res, 1, ScalarInteger(status));
    SET_VECTOR_ELT(res, 2, mkString(phase != NULL ? phase : "unknown"));
    SET_VECTOR_ELT(res, 3, output_fastq != NULL ? mkString(output_fastq) : ScalarString(NA_STRING));
    SET_VECTOR_ELT(res, 4, log_file != NULL ? mkString(log_file) : ScalarString(NA_STRING));
    SET_VECTOR_ELT(res, 5, mkString(stderr_text != NULL ? stderr_text : ""));

    setAttrib(res, R_NamesSymbol, names);

    UNPROTECT(2);
    return res;
}

SEXP pandaseq_merge_tokens(SEXP args_in)
{
    int argc, i;
    char **argv;

    PandaArgsFastq data;
    PandaAssembler assembler;
    PandaMux mux;
    PandaOutputSeq output = NULL;
    void *output_data = NULL;
    PandaDestroy output_destroy = NULL;
    int threads = 1;

    stderr_capture_t cap;
    char *parse_stderr = NULL;
    char *run_stderr = NULL;
    char *all_stderr = NULL;

    bool parse_ok;
    bool run_ok;

    const char *output_path;
    const char *log_path;

    if (!Rf_isString(args_in)) {
        Rf_error("args must be a character vector");
    }

    argc = LENGTH(args_in);
    if (argc < 1) {
        Rf_error("args must contain at least one token");
    }

    argv = (char **) R_alloc((size_t) argc, sizeof(char *));
    for (i = 0; i < argc; i++) {
        argv[i] = (char *) CHAR(STRING_ELT(args_in, i));
    }

    output_path = find_option_value(argc, argv, "-w");
    log_path = find_option_value(argc, argv, "-g");

    data = panda_args_fastq_new();
    if (data == NULL) {
        Rf_error("Failed to allocate PandaArgsFastq");
    }

    if (stderr_capture_begin(&cap) != 0) {
        panda_args_fastq_free(data);
        Rf_error("Failed to capture stderr for PANDAseq parser");
    }

    parse_ok = panda_parse_args(
        argv, argc,
        panda_stdargs, panda_stdargs_length,
        panda_args_fastq_args, panda_args_fastq_args_length,
        (PandaTweakGeneral) panda_args_fastq_tweak,
        (PandaOpener) panda_args_fastq_opener,
        (PandaSetup) panda_args_fastq_setup,
        data, &assembler, &mux, &threads,
        &output, &output_data, &output_destroy
    );

    parse_stderr = stderr_capture_end(&cap);

    if (!parse_ok) {
        panda_args_fastq_free(data);
        return make_result(
            0,
            1,
            "parse",
            output_path,
            log_path,
            parse_stderr
        );
    }

    if (stderr_capture_begin(&cap) != 0) {
        panda_args_fastq_free(data);
        free(parse_stderr);
        Rf_error("Failed to capture stderr for PANDAseq run");
    }

    run_ok = panda_run_pool(threads, assembler, mux, output, output_data, output_destroy);

    run_stderr = stderr_capture_end(&cap);

    panda_args_fastq_free(data);

    all_stderr = join_messages(parse_stderr, run_stderr);
    free(parse_stderr);
    free(run_stderr);

    return make_result(
        run_ok ? 1 : 0,
        run_ok ? 0 : 2,
        run_ok ? "ok" : "run",
        output_path,
        log_path,
        all_stderr
    );
}