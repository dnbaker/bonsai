.PHONY=all tests clean bonsai update

all:
	$(MAKE) -C bonsai

clean:
	$(MAKE) clean -C bonsai

unit:
	$(MAKE) update unit -C bonsai

zunit:
	$(MAKE) update zunit -C bonsai

bonsai:
	$(MAKE) update bonsai -C bonsai

bonsai_d:
	$(MAKE) update bonsai_d -C bonsai

distcmp:
	$(MAKE) update distcmp -C bonsai

update:
	$(MAKE) update -C bonsai
