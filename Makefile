.PHONY=all tests clean bonsai update

all:
	$(MAKE) -C bonsai bonsai_z

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

bonsai_z:
	$(MAKE) update bonsai_z -C bonsai

distcmp_z:
	$(MAKE) update distcmp_z -C bonsai

distcmp:
	$(MAKE) update distcmp -C bonsai

update:
	$(MAKE) update -C bonsai
