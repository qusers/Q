#ifndef __Q_INPUT_H__
#define __Q_INPUT_H__

enum q_input_mode_t {
    Q_INPUT_CSV = 0,
    Q_INPUT_NATIVE = 1
};

extern q_input_mode_t q_input_mode;
extern bool q_native_fresh_start;

void q_prepare_input(const char *input_file, const char *workdir);

#endif /* __Q_INPUT_H__ */
