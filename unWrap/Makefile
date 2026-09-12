CC =		gcc
MACHTYPE =	$(shell uname -m)
OSTYPE =	$(shell uname -s)
OBJDIR =	$(MACHTYPE)-$(OSTYPE)

$(shell mkdir -p $(OBJDIR))

CFLAGS =	$(FLAGS) -c -I$(INCLUDEPATH) -MMD -MP

SRCS =	unwrap.c unwrapPhase.c phaseImage.c residueImage.c branchCuts.c \
	addPhaseRamp.c bcut.c labelRegions.c interpolatePhase.c

OBJS =	$(SRCS:%.c=$(OBJDIR)/%.o)

# all: MUST come before the -include. Each .d declares a rule, and the first
# rule make sees becomes the default goal, so with the include first a bare
# "make" (which is what the parent Makefile issues) builds only the first .d's
# object and silently skips every other source in this directory.
all:	$(OBJS)

-include $(OBJS:.o=.d)

.KEEP_STATE:

$(OBJDIR)/%.o: %.c
	$(CC) $(CFLAGS) $< -o $@
