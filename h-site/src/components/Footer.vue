<template>
  <v-card>
    <v-footer padless relative class="font-weight-small">
      <v-col class="text-center" cols="6">
        Created by:
        <p>Halvor Nafstad</p>
      </v-col>
      <v-col class="text-center">
        Contact:
        <p>mrnafstad@gmail.com</p>
      </v-col>
    </v-footer>
  </v-card>
</template>

<script>
export default {
  name: "Footer"
};
</script>

<style scoped>
.v-application .text-center {
  font-size: 75%;
}
p {
  margin-bottom: 0;
}
</style>