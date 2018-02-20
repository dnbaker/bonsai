.PHONY=all tests clean

all:
	$(MAKE) -C bonsai

clean:
	$(MAKE) clean -C bonsai

unit:
	$(MAKE) unit -C bonsai

zunit:
	$(MAKE) zunit -C bonsai

bonsai:
	$(MAKE) bonsai -C bonsai

bonsai_z:
	$(MAKE) bonsai_z -C bonsai

distcmp_z:
	$(MAKE) distcmp_z -C bonsai

distcmp:
	$(MAKE) distcmp -C bonsai
