# pil - Textbook Physics Illustration Library
# core (libm only) + Cairo renderer, optional SDL2 window
# make [test] [install PREFIX=/usr/local] [clean]

CC      ?= cc
CFLAGS  ?= -O2 -g
PREFIX  ?= /usr/local
BUILD   := build

override CFLAGS += -std=gnu99 -Wall -Wextra -Iinclude -Isrc $(shell pkg-config --cflags cairo)
override LIBS := $(shell pkg-config --libs cairo) -lm
SRCS    := src/world.c src/elements.c src/solver.c src/draw.c

HAVE_SDL := $(shell pkg-config --exists sdl2 && echo 1)
ifeq ($(HAVE_SDL),1)
  override CFLAGS += -Dpil_HAVE_SDL $(shell pkg-config --cflags sdl2)
  override LIBS := $(shell pkg-config --libs sdl2) $(LIBS)
  SRCS   += src/app_sdl.c
endif

OBJS     := $(SRCS:src/%.c=$(BUILD)/%.o)
EXAMPLES := $(patsubst examples/%.c,$(BUILD)/%,$(filter-out examples/example_main.h,$(wildcard examples/*.c)))

all: $(BUILD)/libpil.a $(EXAMPLES) $(BUILD)/test_core
	@echo "built: libpil.a, $(words $(EXAMPLES)) examples, tests   (SDL2 window: $(if $(HAVE_SDL),yes,no))"

$(BUILD)/%.o: src/%.c $(wildcard src/*.h include/pil/*.h) | $(BUILD)
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILD)/libpil.a: $(OBJS)
	ar rcs $@ $^

$(BUILD)/%: examples/%.c examples/example_main.h $(BUILD)/libpil.a
	$(CC) $(CFLAGS) -Iexamples $< $(BUILD)/libpil.a $(LIBS) -o $@

$(BUILD)/test_core: tests/test_core.c $(BUILD)/libpil.a
	$(CC) $(CFLAGS) $< $(BUILD)/libpil.a $(LIBS) -o $@

$(BUILD):
	mkdir -p $(BUILD)

test: $(BUILD)/test_core
	./$(BUILD)/test_core

install: $(BUILD)/libpil.a
	install -d $(PREFIX)/lib $(PREFIX)/include/pil $(PREFIX)/lib/pkgconfig
	install -m 644 $(BUILD)/libpil.a $(PREFIX)/lib/
	install -m 644 include/pil/*.h $(PREFIX)/include/pil/
	@printf 'prefix=$(PREFIX)\nName: pil\nDescription: textbook physics\nVersion: 0.1\nRequires: cairo\nCflags: -I$${prefix}/include\nLibs: -L$${prefix}/lib -lpil -lm $(if $(HAVE_SDL),-lSDL2)\n' > $(PREFIX)/lib/pkgconfig/pil.pc

clean:
	rm -rf $(BUILD)

.PHONY: all test install clean
