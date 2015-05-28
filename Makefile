########################################################################
# Main binary file
TARGET   := pe

# Directories
BUILD    := build
INCLUDES := include
SOURCES  := src
TESTS    := tests
LIBDIRS  :=

# Compilation flags
CFLAGS  := -g -std=c99 -Wall $(INCLUDE)
NVFLAGS := -g $(INCLUDE)
LDFLAGS := -g

# Additional libraries
LIBS := -lm -lblas -llapacke
########################################################################

.SUFFIXES:

ifneq ($(BUILD),$(notdir $(CURDIR)))

export CC := gcc
#export NV := nvcc
export LD := $(CC)

export OUTPUT := $(CURDIR)/$(TARGET)

export VPATH := $(CURDIR)/$(subst /,,$(dir $(ICON))) \
                $(foreach dir,$(SOURCES),$(CURDIR)/$(dir)) \
                $(foreach dir,$(TESTS),$(CURDIR)/$(dir))

export DEPSDIR := $(CURDIR)/$(BUILD)

CFILES  := $(foreach dir,$(SOURCES),$(notdir $(wildcard $(dir)/*.c)))
#CUFILES := $(foreach dir,$(SOURCES),$(notdir $(wildcard $(dir)/*.cu)))
export OFILES   := $(CFILES:.c=.o) #$(CUFILES:.cu=.o)

TSTCFILES := $(foreach dir,$(TESTS),$(notdir $(wildcard $(dir)/*.c)))
export OTSTFILES := $(TSTCFILES:.c=.o)
export TSTOUTPUT := $(TSTCFILES:.c=.tst)

export INCLUDE  := $(foreach dir,$(INCLUDES),-I $(CURDIR)/$(dir)) \
                   $(foreach dir,$(LIBDIRS),-I$(dir)/include) \
                   -I$(CURDIR)/$(BUILD)
export LIBPATHS := $(foreach dir,$(LIBDIRS),-L$(dir)/lib)

.PHONY: $(BUILD) clean

$(BUILD):
	@mkdir -p $@
	@$(MAKE) --no-print-directory -C $(BUILD) -f $(CURDIR)/Makefile

clean:
	@echo [RM] $(BUILD) $(OUTPUT)
	@rm -fr $(BUILD) $(OUTPUT)

else

all: $(OUTPUT) $(OTSTFILES) $(TSTOUTPUT)

$(OUTPUT): $(OFILES)
	@echo [LD] $(notdir $@)
	@$(LD) $(LDFLAGS) $(OFILES) $(LIBPATHS) $(LIBS) -o $@

%.tst: %.o
	@echo [LD] $(notdir $@)
	@$(LD) $(LDFLAGS) $< $(filter-out main.o,$(OFILES)) \
	       $(LIBPATHS) $(LIBS) -o $@

%.o: %.c
	@echo [CC] $(notdir $<)
	@$(CC) -MMD -MP -MF $(DEPSDIR)/$*.d $(CFLAGS) -c $< -o $@

#%.o: %.cu
#	@echo [NV] $(notdir $<)
#	@$(NV) $(NVFLAGS) -c $< -o $@

-include $(DEPSDIR)/*.d

endif
