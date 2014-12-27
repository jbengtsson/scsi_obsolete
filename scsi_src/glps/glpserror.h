#ifndef GLPS_ERRNO_H
#define GLPS_ERRNO_H

#ifdef __cplusplus
extern "C" {
#endif

enum {
    GLPS_SUCCESS  = 0,
    GLPS_FAILURE  = -1,
    GLPS_CONTINUE = -2,
    GLPS_EEOF     = 1,  /* end of file */
    GLPS_ENOLAT   = 2,  /* lattice file not found */ 
    GLPS_ESYNTAX  = 3,  /* syntax error */
    GLPS_EBLNODEF = 4,  /* beamline is not defined */
    GLPS_ENOVAR,        /* variable not defined */
    GLPS_ECALLABLE,     /* variable is not callable */
    GLPS_ENOFUNC,       /* function is not defined */
    GLPS_EFUNCPAR,      /* function parameters error */
    GLPS_ENOELEMT,      /* element not defined */
    GLPS_ENOPRPT,       /* property is not defined in element */ 
    GLPS_ENOACT,        /* action is not defined */
    GLPS_ENOBL,         /* no beamline defined in lattice file */
    GLPS_ENOTREAL,      /* not a real number */
    GLPS_ENOTVEC,       /* not a vector */
    GLPS_ENOTSTR,       /* not a string */
    GLPS_EGRAMMAR,      /* grammar error */
    GLPS_EBOUNDARY,     /* out of boundary */
    GLPS_EASSIGN,       /* incompatible assignment */
    GLPS_EFLAT,         /* can not write flat file */
    GLPS_EXML,          /* can not write xml file as output */
    GLPS_EZERODIV,      /* divided by zero */
    GLPS_EFOPEN,        /* open file error */
    GLPS_ECHAR,         /* unknown character */
    GLPS_EERR     =256  /* general error */
};

 #define GLPS_WARN(msg, args ...)               \
    fprintf(stderr, msg, args)

void glps_stream_printf (const char *label, const char *file,
                         int line, const char *reason);

const char * glps_strerror (const int gsl_errno);


/* GLPS_ERROR: call the error handler, and return the error code */

#define GLPS_ERROR(reason, glps_errno) \
    do {                                                     \
        glps_error (reason, __FILE__, __LINE__, glps_errno) ;   \
        return glps_errno ;                                     \
    } while (0)

/* GLPS_ERROR_VAL: call the error handler, and return the given value */

#define GLPS_ERROR_VAL(reason, file, line, lat, glps_errno, value)      \
    do {                                                     \
        glps_lat_error (reason, file, line, lat, glps_errno) ; \
        return value ;                                          \
    } while (0)

#define GLPS_LAT_ERROR(reason, file, line, lat, glps_errno)     \
    do { \
        glps_lat_error(reason, file, line, lat, glps_errno) ;    \
        return glps_errno;                                       \
    } while (0)
    


#define DBG_MSG(msg, args ...) \
 fprintf(stderr, "[%s:%d:%s]: " msg, __FILE__, __LINE__, __func__, ##args)

#ifndef NDEBUG
#define GLPS_ERROR_MSG(code, msg, args ...)                 \
 DBG_MSG(msg, ##args);                                                    \
 glps_err_code = code;                                                  \
 sprintf(glps_err_msg, "=glps= [ERROR]: %s, line %d-%d, column %d-%d\n" \
        "=glps=   " msg, glps_current_latfile,                          \
         yylloc.first_line, yylloc.last_line,                           \
         yylloc.first_column, yylloc.last_column, ##args);              \
 glps_error(code, glps_err_msg)
#else
#define GLPS_ERROR_MSG(code, msg, args ...)                 \
 glps_err_code = code;                                                  \
 sprintf(glps_err_msg, "=glps= [ERROR]: %s, line %d-%d, column %d-%d\n" \
         "=glps=   " msg, glps_current_latfile,                         \
         yylloc.first_line, yylloc.last_line,                           \
         yylloc.first_column, yylloc.last_column, ##args);              \
 glps_error(code, glps_err_msg) 
#endif

#ifdef __cplusplus
}
#endif
#endif
