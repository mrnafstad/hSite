<template>
  <v-container class="text-center">
    This is an app for tracking bikes in Oslo
    <Map />
  </v-container>
</template>

<script>
// @ is an alias to /src
import Map from '../components/Map.vue'
export default {
  name: 'Home',
  components: {
    Map
  }
}
</script>
