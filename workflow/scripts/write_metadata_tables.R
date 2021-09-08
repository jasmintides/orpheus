#!/usr/bin/env Rscript

seed_no = 1
char_len = 5

create_guid = function(seed_no = 1, char_len = 14) {
  pool = c(letters, LETTERS, 0:9)
  guid = paste0(sample(pool, char_len, replace = T), collapse = "")
  return (guid)
}

create_guid()

fact_table = data.frame(
  agiRxiv_id = create_guid(),
  orpheus_id = orpheus_id
)

agiRxiv_table = data.frame()
